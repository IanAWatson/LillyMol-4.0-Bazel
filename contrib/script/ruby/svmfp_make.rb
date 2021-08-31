# Build svmfp models

# frozen_string_literal: true

require 'date'
require 'fileutils'
require 'google/protobuf'

require_relative 'lib/iwcmdline'
require_relative 'lib/svmfp_model_pb'

def usage(retcod)
  $stderr << "Builds an svmfp model from smiles and activity\n"
  $stderr << " -mdir <dir>   model directory to create\n"
  $stderr << " -A <fname>    file containing activity data\n"
  $stderr << " -C            classification model\n"
  $stderr << " -gfp ... -gfp fingerprint specification (to gfp_make)\n"
  $stderr << " -p <support>  support level for bit inclusion\n"
  $stderr << " -v            verbose output\n";

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
  cmd = "class_label_translation -C #{mdir}/class_label_translation #{activity_file} > #{train_activity}"
  execute_cmd(cmd, verbose, [train_activity])
end

def get_response_name(activity_file)
  input = File.open(activity_file, 'r')
  header = input.readline
  input.close
  header.split[1]
end

cmdline = IWCmdline.new('-v-mdir=s-A=sfile-C-gfp=close-svml=close-p=ipos-gfp_make=xfile' + 
                        '-svm_learn=xfile-gfp_to_svm_lite=xfile')
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

if ARGV.empty?
  $stderr << "Insufficient arguments\n"
  usage(1)
end

fingerprints = if cmdline.option_present('gfp')
                 cmdline.value('gfp')
               else
                 '-EC3:ACHRY'
               end

smiles = ARGV[0]
activity_file = cmdline.value('A')

train_smi = "#{mdir}/train.smi"
train_gfp = "#{mdir}/train.gfp"
train_activity = "#{mdir}/train.activity"

cmd = "#{gfp_make} #{fingerprints} #{smiles} > #{train_gfp}"
execute_cmd(cmd, verbose, [train_gfp])

FileUtils.cp(smiles, train_smi)
if cmdline.option_present('C')
  perform_class_label_translation(activity_file, mdir, train_activity, verbose)
  svm_learn_options = svm_learn_options + ' -z c'
else
  FileUtils.cp(activity_file, train_activity)
  svm_learn_options = svm_learn_options + ' -z r'
end

bit_xref = "#{mdir}/bit_xref"
cmd = "#{gfp_to_svm_lite} -C #{bit_xref} -A #{train_activity}"
cmd << ' -p ' << cmdline.value('p') if cmdline.option_present('p')

train_svml = "#{mdir}/train.svml"
cmd = "#{cmd} #{train_gfp} > #{train_svml}"
execute_cmd(cmd, verbose, [train_svml])

model_file = "#{mdir}/train.model"
cmd = "#{svm_learn} #{svm_learn_options} #{train_svml} #{model_file}"
execute_cmd(cmd, verbose, [model_file])

support_vectors = "#{mdir}/support_vectors.gfp"
cmd = "svm_model_support_vectors.sh -o #{support_vectors} #{model_file} #{train_gfp}"
execute_cmd(cmd, verbose, [support_vectors])

# Now create the model proto

model = SvmfpModel::SvmfpModel.new
model.threshold_b = get_threshold_b(model_file)
model.bit_xref = 'bit_xref'
model.fingerprints = fingerprints
model.support_vectors = 'support_vectors.gfp'
model.date_built = Time.now.to_s
model.response_name = get_response_name(train_activity)
model.class_label_translation = 'class_label_translation' if cmdline.option_present('C')

model_description_file = "#{mdir}/model.dat"
File.write(model_description_file, SvmfpModel::SvmfpModel.encode(model))
model_description_file = "#{mdir}/model.json"
File.write(model_description_file, SvmfpModel::SvmfpModel.encode_json(model))
