#!/usr/bin/env bash

VENV_PATH='etomica-core/src/main/resources/virtualenv'

python -m venv ${VENV_PATH}

source "${VENV_PATH}/bin/activate"

pip install -r config/requirements.txt
