# Build a descriptor model.
# Takes a smiles file and an activity file, together with a list
# of commands that can generate descriptors, iwdescr.sh by default
# frozen_string_literal: true

require 'date'
require 'fileutils'
require 'tempfile'
require 'google/protobuf'

require_relative 'lib/iwcmdline'
require_relative 'lib/gfp_model_pb'

def usage
  $stderr << "Build a molecular descriptor model\n"
  $stderr << " -A <fname>   activity file data (mandatory)\n"
  $stderr << " -mdir <dir>  create model in <mdir>\n"
  $stderr << " -xgboost ... -xgboost  build an xgboost model\n"
  $stderr << " -xgboost_config <fname> xgboost config file\n"
  $stderr << " -desc <script>  specify one or more tools that generate descriptors\n"
  $stderr << " -c            classification model\n"
  exit
end

def execute_cmd(cmd, verbose, files_created)
  $stderr << "Executing '#{cmd}'\n" if verbose
  system(cmd)
  files_created.each do |fname|
    unless File.size?(fname)
      $stderr << "#{cmd} did not create #{fname}\n"
      exit
    end
  end

  true
end

# Compute one or more descriptor sets for `smiles`.
def compute_descriptors(verbose, smiles, scripts, destination)
  if scripts.length == 1
    cmd = "#{scripts.first} #{smiles} > #{destination}"
    execute_cmd(cmd, verbose, [destination])
    return true
  end

  unlink_now = false
  # Multiple scripts, we just execute each and join the resulting files
  tempfiles = []
  scripts.each do |script|
    tmp = Tempfile.new('compute_descriptors')
    tmp.close(unlink_now)
    cmd = "#{script} #{smiles} > #{tmp.path}"
    execute_cmd(cmd, verbose, [tmp.path])
    tempfiles.push(tmp.path)
  end

  cmd = "concat_files #{tempfiles.join(' ')} > #{destination}"
  execute_cmd(cmd, verbose, [destination])
  true
end

# Convert `descriptor_file` to libsvm form as `destination_fname`.
# `activity_fname` is a descriptor file holding the activity data.
# We return the response name.
def convert_to_libsvm(verbose, descriptor_file, activity_fname, destination_fname)
  cmd = "descriptor_file_to_svml -A #{activity_fname} #{descriptor_file} > #{destination_fname}"
  return nil unless execute_cmd(cmd, verbose, [destination_fname])

  # Fetch first line from file.
  header = File.open(activity_fname, "r").first
  return header.split[0]

# This works, but c++ based tool is better.
#  activity = {}
#  response_name = ''
#  File.readlines(activity_fname).each do |line|
#    if response_name.empty?
#      response_name = line.chomp.split[1]
#      next
#    end
#    f = line.chomp.split
#    activity[f[0]] = f[1]
#  end
#
#  $stderr << "Read #{activity.size} activity values from #{activity_fname}\n" if verbose
#  File.open(destination_fname, 'w') do |destination|
#    first_line = true
#    File.readlines(descriptor_file).each do |line|
#      if first_line
#        first_line = false
#        next
#      end
#      f = line.chomp.split
#      unless activity.key?(f[0])
#        $stderr << "No activity for #{f[0]}\n"
#        return false
#      end
#      destination << activity[f[0]]
#      # $stderr << "id #{f[0]} activity #{activity[f[0]]}\n"
#
#      f[1..].each_with_index do |token, ndx|
#        destination << " #{ndx + 1}:#{token}"
#      end
#      destination << "\n"
#    end
#  end
#
# response_name
end

# A model file has been created in `mdir`. return the name.
# Note that this will silently fail of there are multiple model
# files in the directory and the one we want is not the first found.
def get_xgboost_model_file(mdir)
  maybe = Dir.children(mdir).filter { |fname| /^\d+\.model$/.match(fname) }

  if maybe.empty?
    $stderr << "No model file in #{mdir}\n"
    return nil
  end

  return maybe.first
end

# The metadata attribute is common among all model types.
def populate_metadata(model, response_name, classification)
  model.metadata = GfpModel::ModelMetadata.new
  model.metadata.date_built = Time.now.to_s
  model.metadata.response_name = response_name
  if classification
    model.metadata.class_label_translation = 'class_label_translation.dat'
  else
    model.metadata.response_scaling = 'response_scaling.dat'
  end
end

def main
  cl = IWCmdline.new('-v-A=sfile-mdir=s-xgboost=close-xgboost_config=sfile-desc=s-C')

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage
  end

  verbose = cl.option_present('v')

  unless cl.option_present('A')
    $stderr << "Must specity activity file via the -A option\n"
    usage
  end

  unless cl.option_present('mdir')
    $stderr << "must specify model directory via the -mdir option\n"
    usage
  end

  if ARGV.size != 1
    $stderr << "Must specify smiles file as an argument\n"
    usage
  end

  unless cl.option_present('xgboost_config')
    $stderr << "Must specify xgboost config file via the -xgboost_config option\n"
    usage
  end

  xgboost_config = cl.value('xgboost_config')

  smiles = ARGV[0]

  mdir = cl.value('mdir')

  begin
    FileUtils.mkdir_p(mdir)
  rescue => e # rubocop:disable Style/RescueStandardError
    $stderr << "Did not create model directory #{mdir} #{e.message}'\n"
    exit(1)
  end

  scripts = cl.values('desc')
  scripts << 'iwdescr -g all -l -O all' if scripts.empty?

  descriptors = Tempfile.new('descriptor_model')
  dpath = descriptors.path
  unlink_now = false
  descriptors.close(unlink_now)
  compute_descriptors(verbose, smiles, scripts, dpath)
  unless File.size?(dpath)
    $stderr << "Descriptor file #{dpath} not generated\n"
    exit 1
  end

  libsvm = "#{mdir}/train.libsvm"
  response_name = convert_to_libsvm(verbose, dpath, cl.value('A'), libsvm)

  cmd = String.new
  cmd << 'xgboost'
  cmd << ' ' << xgboost_config
  cmd << ' ' << cl.value('xgboost') if cl.option_present('xgboost')
  cmd << if cl.option_present('C')
            ' objective=binary:logistic'
          else
            ' objective=reg:squarederror'
          end
  cmd << " model_dir=#{mdir}"
  cmd << " data=#{libsvm}?format=libsvm"
  execute_cmd(cmd, verbose, [])

  model = GfpModel::DescriptorModel.new
  populate_metadata(model, response_name, cl.option_present('C'))
  scripts.each do |script|
    model.descriptor_generator << script
  end

  model.model_file = get_xgboost_model_file(mdir)

  File.write("#{mdir}/model.dat", GfpModel::DescriptorModel.encode(model))
  File.write("#{mdir}/model.json", GfpModel::DescriptorModel.encode_json(model))
  FileUtils.cp(xgboost_config, "#{mdir}/xgboost.config")
  FileUtils.cp(smiles, mdir)
  FileUtils.cp(cl.value('A'), "#{mdir}/train.activity")
end

main
