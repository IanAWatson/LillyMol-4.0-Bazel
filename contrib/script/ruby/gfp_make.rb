#!/usr/bin/env ruby

# frozen_string_literal: true

# Fingerprint generator

require_relative 'lib/fp_config_files'
require_relative 'lib/group_args'
require_relative 'lib/mkmk2'

def usage(config_dirs, retcode)
  fps = config_fingerprints(config_dirs, 0)
  $stderr << "Generates fingerprints\n"
  $stderr << "The following fingerprints are recognised\n"
  $stderr << fps.keys.sort.map { |s| "-#{s}" }.join(' ') << "\n"
  $stderr << "Many fingerprint generators can be qualified with a radius or distance, -EC3, -AP9, -MK2\n"
  $stderr << "To pass options to a fingerprint generator an opening and closing option are needed\n"
  $stderr << "   gfp_make -oFOO -a 1 -B foo -cFOO ...\n"
  $stderr << " will pass '-a 1 -B foo' to the FOO fingerprint generator\n"
  $stderr << " -v               verbose output\n"
  exit(retcode)
end

config_dirs = ["#{__dir__}/FP"]

if ARGV.empty?
  $stderr << "No arguments\n"
  usage(config_dirs, 1)
end

# First scan arguments to extract options applicable this file.

verbose = 0
is_filter = false
remaining_args = []
argptr = 0
while argptr < ARGV.size
  arg = ARGV[argptr]
  argptr += 1
  if arg == '-v'
    verbose += 1
    next
  end
  if /^--?DIR$/.match(arg)
    config_dirs.push!(ARGV[argptr])
    argptr += 1
    next
  end
  m = /-DIR=(\S+)/.match(arg)
  if m
    config_dirs.push!(m[1])
    next
  end
  if arg == '-f'
    is_filter = true
    next
  end
  remaining_args.push(arg)
end

$stderr << "Config dirs #{config_dirs}\n" if verbose.positive?
config_dirs.each do |dir|
  unless Dir.exist?(dir)
    $stderr << "Config directory #{dir} not found\n"
    exit(1)
  end
end

if remaining_args.empty?
  $stderr << "No fingerprint arguments or inputs\n"
  usage(config_dirs, 1)
end

# By convention, input file(s) must follow arguments.
files = []
while remaining_args.size.positive?
  arg = remaining_args.last
  if arg == '-'
    remaining_args.pop
    files.push(arg)
    next
  end

  break if arg.start_with?('-')

  files.push(remaining_args.pop)
end
files = files.reverse

fp_args = GfpMakeSupport.group_args(remaining_args)
unless fp_args
  fp_args = [OptionValue.new('IW', nil), OptionValue.new('MK', '-J LEVEL2=MK2')]
  $stderr << "Cannot parse fingerprint arguments\n"
  $stderr << "default fp #{fp_args}\n"
  # exit(1)
end

# There is one piece of Magic needed. If -MK and -MK2 are specified,
# combine them.

fp_args = consoliate_mkmk2(fp_args)

# Detect any duplicates

unique_fps = GfpMakeSupport.all_options_unique(fp_args)
unless unique_fps
  $stderr << 'Duplicate fingerprints detected\n'
  exit(1)
end

# Keep track of the recognized fingerprints.
# A mapping from fingerprint name to an object that
# can process that kind of fingerprint.
fps = config_fingerprints(config_dirs, verbose)

# For each fp_option, a mapping to an object that knows how to
# generate that command line component.
fp_option_to_known_fp = {}

fp_args.each do |fp_option|
  matched = 0
  opt = fp_option.option
  fps.each_value do |fp_generator|
    if fp_generator.match?(opt)
      fp_option_to_known_fp[opt] = fp_generator
      matched += 1
    end
  end
  if matched.zero?
    $stderr << "Unrecognized fingerprint #{opt}\n"
    exit(1)
  elsif matched > 1
    $stderr << "#{matched} known fingerprints match #{opt}\n"
    exit(1)
  elsif verbose.positive?
    $stderr << "Recognized #{opt} as #{fp_option_to_known_fp[opt]}\n"
  end
end

# Now that each of unique_fps is associated with an object, the
# command line can be built.

cmdline = []
first_token = true
first_token = false if is_filter
files = files.join(' ')

fp_args.each do |fp_option|
  cmdline.push(fp_option_to_known_fp[fp_option.option].expand(fp_option.option,
                                                              first_in_pipeline: first_token,
                                                              extra_qualifiers: fp_option.value))
  if first_token
    first_token = false
    cmdline[-1] << " #{files}"
  else
    cmdline[-1] << ' -'
  end
end
$stderr << 'Executing ' << cmdline.join('|') << "\n" if verbose.positive?

system(cmdline.join('|'))
