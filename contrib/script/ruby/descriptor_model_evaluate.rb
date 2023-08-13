# Evaluate an xgboost descriptor model built with descriptor_model_build

require 'tempfile'

require_relative 'lib/iwcmdline.rb'
require_relative 'lib/gfp_model_pb'

def usage
  $stderr << "Score a molecular descriptor model built with descriptor_model_make\n"
  $stderr << " -mdir <dir>  model directory\n"
  $stderr << " -xgboost ... -xgboost  build an xgboost model\n"
end

def execute_cmd(cmd, verbose, files_created)
  $stderr << "Executing '#{cmd}\n" if verbose
  system(cmd)
  files_created.each do |fname|
    unless File.size?(fname)
      $stderr << "#{cmd} did not create #{fname}\n"
      exit
    end
  end

  true
end

# Return the proto contents of "#{mdir}/model.dat"
def get_model_proto(fname)
  raise "#{fname} not present" unless File.size?(fname)
  GfpModel::DescriptorModel.decode(File.read(fname))
end

# Compute one or more descriptor sets for `smiles`.
def compute_descriptors(verbose, smiles, scripts, destination)
  if scripts.length == 1
    cmd = "#{scripts.first} #{smiles} > #{destination}"
    execute_cmd(cmd, verbose, [destination])
    return true
  end

  # Multiple scripts, we just execute each and join the resulting files
  tempfiles = []
  scripts.each do |script|
    tmp = Tempfile.new("compute_descriptors").path
    cmd = "#{script} #{smiles} > #{tmp}"
    execute_cmd(cmd, verbose, [tmp])
    tempfiles.push(tmp)
  end

  cmd = "concat_files #{tempfiles.join(' ')} > #{destination}"
  execute_cmd(cmd, verbose, [destination])
  return true
end

def convert_to_libsvm(verbose, input, destination)
  first_line = true
  File.readlines(input).each do |line|
    if first_line
      first_line = false
      next
    end
    f = line.chomp.split
    destination << "0";
    f[1..].each_with_index do |token, ndx|
      destination << " #{ndx + 1}:#{token}"
    end
    destination << "\n"
  end
end

# return the id's in a smiles file
def get_ids(fname)
  result = []
  File.readlines(fname).each do |line|
    f = line.chomp.split
    result << f[1]
  end
  result
end

def main
  cl = IWCmdline.new('-v-mdir=s-xgboost=close-desc=s')

  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage
  end

  verbose = cl.option_present('v')

  unless cl.option_present('mdir')
    $stderr << "must specify model directory via the -mdir option\n"
    usage
  end

  if ARGV.size != 1
    $stderr << "Must specify smiles file as an argument\n"
    usage
  end

  smiles = ARGV[0]

  mdir = cl.value('mdir')

  unless File.directory?(mdir)
    $stderr << "Model directory #{mdir} not found\n"
    exit 1
  end

  model_proto = get_model_proto("#{mdir}/model.dat")

  xgboost_config = "#{mdir}/xgboost.config"
  unless File.size?(xgboost_config)
    $stderr << "#{xgboost_config} missing or empty\n"
    exit 1
  end

  scripts = model_proto.descriptor_generator

  descriptors = Tempfile.new("descriptor_model_eval")
  dpath = descriptors.path
  descriptors.close(unlink_now=false)
  compute_descriptors(verbose, smiles, scripts, dpath)
  unless File.size?(dpath)
    $stderr << "Descriptor file #{dpath} not generated\n"
    exit 1
  end

  libsvm = Tempfile.new("descriptor_model_eval.libsvm")
  convert_to_libsvm(verbose, dpath, libsvm)
  libsvm.close(unlink_now=false)

  cmd = 'xgboost'
  cmd << " #{xgboost_config}"
  cmd << " test:data=#{libsvm.path}"
  cmd << ' task=pred'
  cmd << " model_in=#{mdir}/#{model_proto.model_file}"
  tmp_pred = Tempfile.new("descriptor_model_pred")
  tmp_pred.close(unlink_now = false)
  cmd << " name_pred=#{tmp_pred.path}"
  execute_cmd(cmd, verbose, [tmp_pred.path])
  ids = get_ids(smiles)
  ndx = 0
  $stdout << "ID #{model_proto.metadata.response_name}\n"
  File.readlines(tmp_pred.path).each do |line|
    $stdout << "#{ids[ndx]} #{line}"
    ndx += 1
  end
end

main
