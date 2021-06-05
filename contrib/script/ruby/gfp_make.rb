#!/usr/bin/env ruby

# frozen_string_literal: true

# Fingerprint generator

require_relative 'lib/fp_config_files'
require_relative 'lib/group_args'
require_relative 'lib/mkmk2'
require 'google/protobuf'

def usage(retcode)
  exit(retcode)
end

if ARGV.empty?
  $stderr << "No arguments\n"
  usage(1)
end

# First scan arguments to extract options for this file.

verbose = 0
is_filter = false
config_dir = "#{__dir__}/FP"
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
    config_dir = ARGV[argptr]
    argptr += 1
    next
  end
  m = /-DIR=(\S+)/.match(arg)
  if m
    config_dir = m[1]
    $stderr << "Config dir set to #{m[1]}\n" if verbose.positive?
    next
  end
  if arg == '-f'
    is_filter = true
    next
  end
  remaining_args.push(arg)
end

$stderr << "Config dir #{config_dir}\n" if verbose.positive?
unless Dir.exist?(config_dir)
  $stderr << "Config directory #{config_dir} not found\n"
  exit(1)
end

if remaining_args.empty?
  $stderr << 'No fingerprint arguments or inputs\n'
  usage(1)
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
  $stderr << 'Cannot parse fingerprint arguments\n'
  exit(1)
end

# There is one piece of Magic needed. If -MK and -MK2 are specified,
# combine them.

fp_args = consoliate_mkmk2(fp_args)

#$stderr << "Files #{files}\n"
#$stderr << "remaining_args #{remaining_args}\n"
#$stderr << "fp_args #{fp_args}\n"

# Detect any duplicates

unique_fps = GfpMakeSupport.all_options_unique(fp_args)
unless unique_fps
  $stderr << 'Duplicate fingerprints detected\n'
  exit(1)
end

#$stderr << "After MK/MK2 consolidation\n"
#$stderr << unique_fps << "\n"

# Keep track of the recognized fingerprints.
# A mapping from fingerprint name to an object that
# can process that kind of fingerprint.
fps = config_fingerprints(config_dir, verbose)

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
$stderr << cmdline.join('|') if verbose.positive?

system(cmdline.join('|'))
