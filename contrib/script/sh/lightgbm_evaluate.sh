#!/bin/sh
dir=$(dirname $0)
ruby "${dir}/../ruby/lightgbm_evaluate.rb" "$@"
