# Build a descriptor model.
# Takes a smiles file and an activity file, together with a list
# of commands that can generate descriptors, iwdescr.sh by default

require 'date'
require 'fileutils'
require 'tempfile'
require 'google/protobuf'

require_relative 'lib/iwcmdline.rb'

def usage
  $stderr << "Build a molecular descriptor model\n"
  $stderr << " -A <fname>   activity file data (mandatory)\n"
  $stderr << " -mdir <dir>  create model in <mdir>\n"
  $stderr << " -xgboost ... -xgboost  build an xgboost model\n"
  $stderr << " -desc <script>  specify one or more tools that generate descriptors\n"
  $stderr << " -c            classification model\n"
  exit
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

def convert_to_libsvm(verbose, descriptors, activity_fname, destination)
  activity = {}
  first_line = true
  File.readlines(activity_fname).each do |line|
    if first_line
      first_line = false
      next
    end
    f = line.chomp.split
    activity[f[0]] = f[1]
  end

  $stderr << "Read #{activity.size} activity values from #{activity_fname}\n" if verbose
  first_line = true
  File.readlines(descriptors).each do |line|
    if first_line
      first_line = false
      next
    end
    f = line.chomp.split
    unless activity.key?(f[0])
      $stderr << "No activity for #{f[0]}\n"
      return false
    end
    destination << activity[f[0]]
    # $stderr << "id #{f[0]} activity #{activity[f[0]]}\n"

    f[1..].each_with_index do |token, ndx|
      destination << " #{ndx + 1}:#{token}"
    end
    destination << "\n"
  end

  true
end

def main
  cl = IWCmdline.new('-v-A=sfile-mdir=s-xgboost=close-desc=s-C')

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

  smiles = ARGV[0]

  mdir = cl.value('mdir')

  begin
    FileUtils.mkdir_p(mdir)
  rescue => e # rubocop:disable Style/RescueStandardError
    $stderr << "Did not create model directory #{mdir} #{e.message}'\n"
    exit(1)
  end

  scripts = cl.values('desc')
  if scripts.empty?
    scripts << 'iwdescr -O all'
  end

  descriptors = Tempfile.new("descriptor_model")
  dpath = descriptors.path
  descriptors.close(unlink_now=false)
  compute_descriptors(verbose, smiles, scripts, dpath)
  unless File.size?(dpath)
    $stderr << "Descriptor file #{dpath} not generated\n"
    exit 1
  end

  libsvm = Tempfile.new("descriptor_model.libsvm")
  convert_to_libsvm(verbose, dpath, cl.value('A'), libsvm)
  libsvm.close(unlink_now=false)

  cmd = "xgboost"
  cmd << ' ' << cl.value('xgboost') if cl.option_present('xgboost')
  if cl.option_present('C')
    cmd << " objective=binary:logistic"
  else
    cmd << " objective=reg:squarederror"
  end
  cmd << " model_dir=#{mdir}"
  cmd << " data=#{libsvm.path}?format=libsvm"
  execute_cmd(cmd, verbose, [])
  #FileUtils.cp

end

main
