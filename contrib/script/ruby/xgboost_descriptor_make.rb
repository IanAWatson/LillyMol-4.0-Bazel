#!/usr/bin/env ruby

# Build an xgboost model using xgboost

require 'fileutils'

require_relative 'lib/iwcmdline'
require_relative 'descriptor_model_pb'

def usage(rc)
  $stderr << "Build an xgboost model using tabular molecular descriptors\n"
  $stderr << " -mdir <dir>      model directory\n"
  $stderr << " -A <fname>       activity file. Header record then 'id activity' pairs\n"
  $stderr << " -C               classification model (no label transalation done here)\n"
  $stderr << " -fmt <csv,svml>  file type passed to xgboost\n"
  $stderr << " -v               verbose output\n"
  exit(rc)
end

def execute_cmd_check(cmd, verbose, check_file)
  $stderr << "Executing #{cmd}\n" if verbose
  unless system(cmd)
    $stderr << "#{cmd} failed\n"
    return false
  end
  check_file.each do |fname|
    unless File.size?(fname)
      $stderr << "#{cmd} did not create #{fname}\n"
      return false
    end
  end
  true
end

def find_model_file(mdir)
  Dir.entries(mdir).each do |fname|
    next if File.directory?(fname)
    return fname if /\.model$/.match(fname)
  end
  $stderr << "No model file in #{mdir}\n"
  None
end

# `header` is the header record in a descriptor file, with the id field dropped.
# Fill the column cross reference data in `proto`
def fill_proto(header, proto)
  header.each_with_index do |token, col|
    p = DescriptorModel::NameToFeature.new

    p.name = token
    p.bit = col
    proto.xref << p
  end
end

# Convert descriptor file `fname` to svml form. The first 'descriptor' in the
# fiie is the activity for that molecule.
def dat_to_svml(fname, proto, svml)
  File.open(svml, 'w') do |writer|
    first_record = true
    IO.foreach(fname).each do |line|
      f = line.chomp.split
      # Discard identifier
      f.shift
      activity = f.shift

      if first_record
        proto.response = activity # First record it is the name of the response
        fill_proto(f, proto)
        first_record = false
        next
      end
      writer << activity
      f.each_with_index do |token, col|
        next if token == '.'

        writer << " #{col}:#{token}"
      end
      writer << "\n"
    end
  end
end

# Write `id_to_activity` to `fname`.
def write_activity(fname, id_to_activity)
  File.open(fname, 'w') do |writer|
    writer << "Id Activity\n"
    id_to_activity.each do |k, v|
      writer << "#{k} #{v}\n"
    end
  end
end

def read_activity_data(fname)
  result = {}
  first_line = true
  IO.foreach(fname).each do |line|
    if first_line
      first_line = false
      next
    end

    f = line.chomp.split
    result[f[0]] = f[1].to_f
  end
  result
end

# Some likely OK hyperparameters for building an xgboost model.
def default_config
  <<~CONFIG
    booster = gbtree
    eta = 0.05
    gamma = 1.0
    min_child_weight = 1
    max_depth = 5
    num_round = 500
  CONFIG
end

def main
  cmdline = IWCmdline.new('-v-A=sfile-mdir=s-xgboost=close-xgboost_config=sfile-C-fmt=s')

  if cmdline.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage(1)
  end

  verbose = cmdline.option_present('v')

  unless cmdline.option_present('mdir')
    $stderr << "Must specify model directory via the -mdir option\n"
    usage(1)
  end

  if cmdline.option_present('A')
    fname = cmdline.value('A')
    id_to_activity = read_activity_data(fname)
    $stderr << "Read #{id_to_activity.size} activity values from #{fname}\n" if verbose
  end

  if ARGV.empty?
    $stderr << "Must specify descriptor file as an argument\n"
    usage(1)
  end

  mdir = cmdline.value('mdir')

  begin
    FileUtils.mkdir_p(mdir)
  rescue => e # rubocop:disable Style/RescueStandardError
    $stderr << "Did not create model directory #{mdir} #{e.message}'\n"
    exit(1)
  end

  config_file = File.join(mdir, 'xgboost.config')
  if cmdline.option_present('xgboost_config')
    FileUtils.cp(cmdline.value('xgboost_config', config_file))
  elif cmdline.option_present('C')
    File.write(config_file, default_config << "\nobjective = binary:logistic\n")
  else
    File.write(config_file, default_config << "\nobjective = reg:squarederror\n")
  end

  train_activity = File.join(mdir, 'train.activity')
  write_activity(train_activity, id_to_activity)

  train_dat = File.join(mdir, 'train.dat')
  cmd = if ARGV.size == 1
          "concat_files #{train_activity} #{ARGV[0]}"
        else
          "concat_files #{train_activity} #{argv.join(' ')}"
        end

  cmd << "> #{train_dat}"

  execute_cmd_check(cmd, verbose, [train_dat])

  train_svml = File.join(mdir, 'train.svml')
  proto = DescriptorModel::Model.new
  dat_to_svml(train_dat, proto, train_svml)

  uri = File.absolute_path(train_svml)

  opts = cmdline.values('xgboost').join(' ')
  cmd = "xgboost #{config_file} model_dir=#{mdir} task=train #{opts} data=#{uri}?format=libsvm"
  execute_cmd_check(cmd, verbose, [])
  model_file = find_model_file(mdir)
  raise 'No model' unless model_file

  proto.model_file = model_file
  File.write("#{mdir}/model.dat", DescriptorModel::Model.encode(proto))
  File.write("#{mdir}/model.json", DescriptorModel::Model.encode_json(proto))

  0
end

main
