#!/bin/bash
dir=$(dirname $0)
exec ruby ${dir}/../ruby/descriptor_model_build.rb "$@"
