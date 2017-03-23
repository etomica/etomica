#!/usr/bin/env bash

VENV_PATH='etomica-core/src/main/resources/virtualenv'

python3 -m venv ${VENV_PATH}

source "${VENV_PATH}/bin/activate"

pip3 install -r config/requirements.txt
