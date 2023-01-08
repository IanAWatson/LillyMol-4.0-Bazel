#!/bin/sh
dir=$(dirname $0)
ruby "${dir}/../ruby/catboost_evaluate.rb" "$@"
