#!/bin/bash
dir=$(dirname $0)

exec ruby ${dir}/../ruby/svmfp_evaluate.rb "$@"
