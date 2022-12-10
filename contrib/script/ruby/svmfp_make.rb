# Build svmfp models
# Subsequently altered to also build lightgbm models.

# frozen_string_literal: true

require 'date'
require 'fileutils'
require 'google/protobuf'

require_relative 'lib/iwcmdline'
require_relative 'lib/gfp_model_pb'

def usage(retcod)
  $stderr << "Builds an svmfp model from smiles and activity\n"
  $stderr << " -mdir <dir>   model directory to create\n"
  $stderr << " -A <fname>    file containing activity data\n"
  $stderr << " -C            classification model\n"
  $stderr << " -gfp ... -gfp fingerprint specification (to gfp_make)\n"
  $stderr << " -p <support>  support level for bit inclusion\n"
  $stderr << " -flatten      flatten sparse fingerprint counts to 1\n"
  $stderr << " -lightgbm ... -lightgbm build a lightgbm model\n"
  $stderr << " -catboost ... -catboost build a catboost model\n"
  $stderr << " -v            verbose output\n"

  exit retcod
end

def execute_cmd(cmd, verbose, files_needed = [])
  $stderr << "Executing '#{cmd}'\n" if verbose
  rc = system(cmd)
  raise "Non zero rc #{rc} from #{cmd}" unless rc

  files_needed.each do |fname|
    next if File.size?(fname)

    raise "#{cmd} did not produce #{fname}"
  end
end

def get_threshold_b(model_file)
  threshold_b_rx = Regexp.new('(\S+) # threshold b, each following line')

  File.foreach(model_file).each do |line|
    m = threshold_b_rx.match(line)
    return m[1].to_f if m
  end
  raise "No threshold_b in #{model_file}"
end

# `activity_file` has been identified as containing classification data.
# generate a class label translation file in `mdir`.
def perform_class_label_translation(activity_file, mdir, train_activity, verbose)
  cmd = "class_label_translation -C #{mdir}/class_label_translation.txt " \
  "-cbin #{mdir}/class_label_translation.dat #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity, "#{mdir}/class_label_translation.dat"])
end

# Only difference is that an extra option is needed for class_label_translation.
def perform_class_label_translation_lightgbm(activity_file, mdir, train_activity, verbose)
  cmd = "class_label_translation -lightgbm -C #{mdir}/class_label_translation.txt " \
  "-cbin #{mdir}/class_label_translation.dat #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity, "#{mdir}/class_label_translation.dat"])
end

def perform_response_scaling(activity_file, mdir, train_smi, train_activity, verbose)
  cmd = "feature_scaling -bin -C #{mdir}/response_scaling -v -subset #{train_smi} -scol 2 #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity])
end

def get_response_name(activity_file)
  input = File.open(activity_file, 'r')
  header = input.readline
  input.close
  header.split[1]
end

cmdline = IWCmdline.new('-v-mdir=s-A=sfile-C-gfp=close-svml=close-p=ipos-flatten-gfp_make=xfile' \
                        '-svm_learn=xfile-gfp_to_svm_lite=xfile-lightgbm=close-lightgbm_config=sfile' \
                        '-catboost=close')
if cmdline.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage(1)
end

verbose = cmdline.option_present('v')

unless cmdline.option_present('A')
  $stderr << "Must specify activity file via the -A option\n"
  usage(1)
end

unless cmdline.option_present('mdir')
  $stderr << "Must specify model directory via the -mdir option\n"
  usage(1)
end

mdir = cmdline.value('mdir')

begin
  FileUtils.mkdir_p(mdir)
rescue => e # rubocop:disable Style/RescueStandardError
  $stderr << "Did not create model directory #{mdir} #{e.message}'\n"
  exit(1)
end

gfp_make = if cmdline.option_present('gfp_make')
             cmdline.value('gfp_make')
           else
             'gfp_make.sh'
           end

gfp_to_svm_lite = if cmdline.option_present('gfp_to_svm_lite')
                    cmdline.value('gfp_to_svm_lite')
                  else
                    'gfp_to_svm_lite'
                  end

svm_learn = if cmdline.option_present('svm_learn')
              cmdline.value('svm_learn')
            else
              'svm_learn'
            end

svm_learn_options = if cmdline.option_present('svml')
                      cmdline.value('svml')
                    else
                      '-t 4 -m 500'
                    end

# nil if not specified.
lightgbm = cmdline.value('lightgbm')
default_lightgbm_config = cmdline.value('lightgbm_config')
catboost = cmdline.value('catboost')

lightgbm = "lightgbm config=#{default_lightgbm_config} #{lightgbm} force_row_wise=true" if lightgbm
catboost = "catboost fit #{catboost} --train-dir #{mdir} --fstr-file fstr.dat " \
           "--use-best-model --min-data-in-leaf=2 " \
           "--model-format CatboostBinary,CPP " if catboost

if lightgbm && ! default_lightgbm_config
  $stderr << "When building a lightgbm model, must specify -lightgbm_config\n"
  usage(1)
end

if ARGV.empty?
  $stderr << "Insufficient arguments\n"
  usage(1)
end

fingerprints = if cmdline.option_present('gfp')
                 cmdline.value('gfp')
               else
                 '-EC3:ACHRY'
               end
flatten_sparse_fingerprints = cmdline.option_present('flatten')

smiles = ARGV[0]
activity_file = cmdline.value('A')

train_smi = "#{mdir}/train.smi"
train_gfp = "#{mdir}/train.gfp"
train_activity = "#{mdir}/train.activity"

cmd = "#{gfp_make} #{fingerprints} #{smiles} > #{train_gfp}"
execute_cmd(cmd, verbose, [train_gfp])

FileUtils.cp(smiles, train_smi)
if cmdline.option_present('C')  # Classification.
  if lightgbm || catboost
    perform_class_label_translation_lightgbm(activity_file, mdir, train_activity, verbose)
    lightgbm = "lightgbm #{lightgbm} objective=binary" if lightgbm
    catboost = "#{catboost} --loss-function Logloss --custom-metric=MCC --auto-class-weights Balanced" if catboost
  else
    perform_class_label_translation(activity_file, mdir, train_activity, verbose)
    svm_learn_options = "#{svm_learn_options} -z c"
  end
else  # Regression
  perform_response_scaling(activity_file, mdir, train_smi, train_activity, verbose)
  svm_learn_options = "#{svm_learn_options} -z r"
  lightgbm = "lightgbm #{lightgbm} objective=regression" if lightgbm
  catboost = "#{catboost} --loss-function RMSE" if catboost
end

bit_xref = "bit"
bit_subset = "bit"

f = if flatten_sparse_fingerprints
      '-f'
    else
      ''
    end

l = if lightgbm || catboost
      '-l'
    else
      ''
    end

cmd = "#{gfp_to_svm_lite} #{f} #{l} -C #{mdir}/#{bit_xref} -A #{train_activity} -S #{mdir}/train "
if cmdline.option_present('p')
  support = cmdline.value('p')
  cmd = "#{cmd} -p #{support}"
end

train_svml = "#{mdir}/train.svml"
cmd = "#{cmd} #{train_gfp}"
execute_cmd(cmd, verbose, [train_svml, "#{mdir}/bit_xref.dat", "#{mdir}/bit_subset.dat"])

if lightgbm
  model_file = "#{mdir}/LightGBM_model.txt"
  cmd = "#{lightgbm} data=#{mdir}/train.svml output_model=#{model_file}"
elsif catboost
  model_file = "#{mdir}/Catboost.model.bin"
  cmd = "#{catboost} --learn-set libsvm://#{mdir}/train.svml --model-file Catboost.model.bin"
else
  model_file = "#{mdir}/train.model"
  cmd = "#{svm_learn} #{svm_learn_options} #{train_svml} #{model_file}"
end
execute_cmd(cmd, verbose, [model_file])

# The metadata attribute is common among all model types.
def populate_metadata(model, fingerprints, response_name, classification, flatten_sparse_fingerprints)
  model.metadata = GfpModel::ModelMetadata.new
  model.metadata.date_built = Time.now.to_s
  model.metadata.fingerprints = fingerprints
  model.metadata.response_name = response_name
  if classification
    model.metadata.class_label_translation = 'class_label_translation.dat'
  else
    model.metadata.response_scaling = 'response_scaling.dat'
  end
  model.metadata.flatten_sparse_fingerprints = flatten_sparse_fingerprints
end

response_name = get_response_name(train_activity)

catboost = "#{catboost} --name=#{response_name}" if catboost

# Write the model proto
if lightgbm
  model = GfpModel::LightGbmModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)
  model.bit_xref = 'bit_xref.dat'
  File.write("#{mdir}/model.dat", GfpModel::LightGbmModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::LightGbmModel.encode_json(model))
elsif catboost
  model = GfpModel::CatboostModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)
  model.bit_xref = 'bit_xref.dat'
  File.write("#{mdir}/model.dat", GfpModel::CatboostModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::CatboostModel.encode_json(model))
else
  support_vectors = "#{mdir}/support_vectors.gfp"
  cmd = "svm_model_support_vectors.sh -o #{support_vectors} #{model_file} #{train_gfp}"
  execute_cmd(cmd, verbose, [support_vectors])

  model = GfpModel::SvmfpModel.new
  populate_metadata(model, fingerprints, response_name, cmdline.option_present('C'), flatten_sparse_fingerprints)

  model.threshold_b = get_threshold_b(model_file)
  model.bit_subset = 'bit_xref.dat'
  model.bit_xref = bit_xref
  model.train_gfp = train_gfp
  model.support_vectors = 'support_vectors.gfp'
  File.write("#{mdir}/model.dat", GfpModel::SvmfpModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::SvmfpModel.encode_json(model))
end
