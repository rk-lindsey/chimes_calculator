#!/bin/bash
set -e

SCRIPT_DIR=$(dirname $0)
SCRIPT_NAME=$(basename $0)

if [ $# -ne 5 ]; then
  echo "Wrong nr. of parameters" >&2
  echo "run_single_test.sh BINARY PARAMETERFILE GEOMETRYFILE CONFIGOPT WORKDIR" >&2
  exit 1
fi

binary=$1
parameterfile=$2
geometryfile=$3
configopt=$4
workdir=$5

if [ ! -d ${workdir} ]; then
  mkdir -p ${workdir}
fi

cd ${workdir}
${binary} \
  ${SCRIPT_DIR}/force_fields/${parameterfile} \
  ${SCRIPT_DIR}/configurations/${geometryfile} \
  ${configopt} >& output

${SCRIPT_DIR}/compare.sh \
  output \
  ${SCRIPT_DIR}/expected_output/${parameterfile}.${geometryfile}.dat
