#!/bin/bash
dir=$(dirname $0)
exec ruby ${dir}/../ruby/build_models.rb "$@"
