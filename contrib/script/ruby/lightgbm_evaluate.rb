# Evaluate a LightGBM model built with svmfp_make

# frozen_string_literal: true

require 'tempfile'
require 'google/protobuf'

require_relative 'lib/iwcmdline'
require_relative 'lib/gfp_model_pb'
require_relative 'lib/class_label_translation_pb'
require_relative 'lib/feature_scaling_pb'

def usage(retcod)
  $stderr << "Evaluate a LightGBM model build with svmfp_make\n"
  $stderr << " -m <dir>     model directory\n"
  $stderr << " -v           verbose output\n"
  exit(retcod)
end

cmdline = IWCmdline.new('-v-mdir=dir-gfp_make=xfile-eval=xfile-eval_opts=close')
if cmdline.unrecognised_options_encountered
  $stderr << "unrecognised_options_encountered\n"
  usage(1)
end

verbose = cmdline.option_present('v')

unless cmdline.option_present('mdir')
  $stderr << "Must specify model directory via the -mdir option\n"
  usage(1)
end

mdir = cmdline.value('mdir')

# A model that will be scored. Attributes are populated by reading
# the proto file in the directory. Note that the proto file is encoded.
class Model
  attr_reader :data
  attr_reader :mdir
  attr_reader :model_data_fname

  def initialize(mdir)
    @mdir = mdir
    @model_data_fname = "#{mdir}/model.dat"
    raise "Missing #{@model_data_fname}" unless File.size?(@model_data_fname)

    @data = GfpModel::LightGbmModel.decode(File.read(@model_data_fname))
    $stderr << @data
  end

  def gfp
    $stderr << "Fingerprints #{@data.metadata.fingerprints}\n"
    @data.metadata.fingerprints
  end

  def is_classification
    @data.metadata.class_label_translation.size.positive?
  end

  def bit_xref
    return "#{@mdir}/#{@data.bit_xref}"
  end
end

# When give a model directory, we decide if that is a single
# model directory, or a directory containing other model directories.

def identify_models(dir)
  model_proto = "#{dir}/model.dat"
  return [dir] if File.size?(model_proto)

  result = []
  Dir.open(dir).each do |fname|
    next if fname.start_with?('.')
    next unless File.directory?(fname)

    model_proto = "#{dir}/model.dat"
    result << dir if File.size?(model_proto)
  end
  result
end

model_dirs = identify_models(mdir)
if model_dirs.empty?
  $stderr << "Invalid model directory #{mdir}\n"
  exit(1)
end

models = []
model_dirs.each do |dir|
  models << Model.new(dir)
end
$stderr << models << "\n"

if ARGV.empty?
  $stderr << "Insufficient arguments\n"
  usage(1)
end

gfp_make = if cmdline.option_present('gfp_make')
             cl.value('gfp_make')
           else
             'gfp_make.sh'
           end

evaluate = if cmdline.option_present('evaluate')
             cl.value('evaluate')
           else
             'lightgbm'
           end

if cmdline.option_present('eval_opts')
  opts = cmdline.value('eval_opts')
  evaluate = "#{evaluate} #{opts}"
end

smiles_files = ARGV.join(' ')

ids = []
ARGF.each do |line|
  f = line.split(' ')
  ids.push(f[1])
end
raise "No smiles" unless ids.size.positive?

# For now, multiple models are not handled.

def execute_cmd(cmd, verbose, files_needed = [])
  $stderr << "Executing '#{cmd}'\n" if verbose
  rc = system(cmd)
  raise "Non zero rc #{rc} from #{cmd}" unless rc

  files_needed.each do |fname|
    next if File.size?(fname)

    raise "#{cmd} did not produce #{fname}"
  end
end


mnames = []
models.each do |model|
  mnames << '-M ' << model.model_data_fname
end
mnames = mnames.join(' ')
$stderr << mnames << "\n"

def paste_files(file1, file2, destination, verbose)
  cmd = "paste -d' ' #{file1} #{file2} > #{destination.path}"
  execute_cmd(cmd, verbose, [destination.path])
  return destination.path
end

def get_class_label_transation(model)
  fname = model.data.metadata.class_label_translation
  input = File.open("#{model.mdir}/#{fname}", 'rb')
  contents = input.read
  return ClassLabelTranslation::ClassLabelTranslation.decode(contents)
end

# Return a hash where the keys are the values
def reverse_hash(hash)
  result = {}
  hash.each do |k, v|
    result[v] = k
  end
  result
end

def write_header(model, classification)
  response = model.data.metadata.response_name

  sep = ' '

  $stdout << "ID#{sep}pred_#{response}"
  $stdout << "#{sep}score" if classification
  $stdout << "\n"
end

def generate_classification_data(model, predictions, ids)
  write_header(model, true)  # true -> classification

  class_label_translation = get_class_label_transation(model)

  sep = ' '

  score_to_label = reverse_hash(class_label_translation.to_numeric)
  n = predictions.size
  n.times do |i|
    score = predictions[i].to_f
    if score <= 0.5
      label = score_to_label[0]
    else
      label = score_to_label[1]
    end
    $stdout << "#{ids[i]}#{sep}#{label}#{sep}#{score}\n"
  end
end

def generate_regression_data(model, predictions, ids)
  write_header(model, false) # false -> not classification

  fname = "#{model.mdir}/#{model.data.metadata.response_scaling}"
  unscaling = FeatureScaling::FeatureScaling.decode(File.read(fname))
  feature_range = unscaling.max - unscaling.min
  sep = ' '
  n = predictions.size
  n.times do |i|
    score = predictions[i].to_f
    score = unscaling.min + score * feature_range
    $stdout << "#{ids[i]}#{sep}#{score.round(3)}\n"
  end
end

models.each do |model|
  tmpsvml = Tempfile.new('lgbm.svml')
  cmd = "#{gfp_make}  #{model.gfp} #{smiles_files} |" \
         "gfp_to_svm_lite -l -X #{model.bit_xref} -S - -O svml - > #{tmpsvml.path}"
  execute_cmd(cmd, verbose, [tmpsvml.path])

  tmplog = Tempfile.new('lgbm.log')
  tmpresult = Tempfile.new('lgbm.txt')
  cmd = "#{evaluate} task=predict input_model=#{mdir}/LightGBM_model.txt " \
        "data=#{tmpsvml.path} predict_disable_shape_check=true " \
        "metric_freq=1000 predict_result=#{tmpresult.path} > #{tmplog.path}"
  execute_cmd(cmd, verbose, [tmpresult.path])
  system("/bin/cat #{tmplog.path} >&2")

  predictions = File.readlines(tmpresult.path)
  raise "Size mismatch #{ids.size} and #{predictions.size}" unless ids.size == predictions.size
  # Generate a combined file
  tmpcombined = Tempfile.new('lgbmb.txt')
  if model.is_classification
    generate_classification_data(model, predictions, ids)
  else
    generate_regression_data(model, predictions, ids)
  end
end
