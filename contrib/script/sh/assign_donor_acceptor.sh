#!/bin/bash
# Assigns donor/acceptor and charges according to Fred Bruns rules.
# Acceptors will get isotope 1, donors isotope 3 and dual mode isotope 2
# This script assumes the queries are located near this script.
# You will need to add '-S output' which will see output.smi being
# created.

script=$(readlink -e $0)
dir=$(dirname ${script})

charges="${dir}/../../data/queries/charges"
hbonds="${dir}/../../data/queries/hbonds"

fileconv -N F:${charges}/queries \
         -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label \
         -g all "$@"
