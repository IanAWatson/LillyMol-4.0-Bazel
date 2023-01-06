#!/bin/bash
dir=$(dirname $0)
exec ruby "${dir}/../ruby/svmfp_make.rb" "$@"
