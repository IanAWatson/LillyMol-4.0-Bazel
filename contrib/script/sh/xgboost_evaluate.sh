#!/bin/sh
dir=$(dirname $0)
ruby "${dir}/../ruby/xgboost_evaluate.rb" "$@"
