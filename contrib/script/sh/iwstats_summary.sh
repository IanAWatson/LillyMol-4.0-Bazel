#!/bin/bash

dir=$(dirname $0)

PYTHONPATH=${dir}/../../.. python "${dir}/../py/iwstats_summary.py" "$@"
