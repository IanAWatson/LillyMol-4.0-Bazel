#!/usr/bin/env ruby
# Evaluate a model built by xgboost_descriptor_make

require 'fileutils'

require_relative 'lib/iwcmdline.rb'
require_relative 'descriptor_model_pb'

def usage(rc)
  $stderr << "Build an xgboost model using tabular molecular descriptors\n"
  exit(1)
end

def main()
  cmdline = IWCmdline.new("-v-mdir=dir")
  unless cmdline.option_present('mdir')
    $stderr << 'Must specify model directory via the -mdir option\n'
    usage(1)
  end
end

main()
