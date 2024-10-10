#!/usr/bin/env bash

SOFTWARE_PATH="/home/flavian/Projects/escape/"

SOFTWARE="$SOFTWARE_PATH/src/main.jl"
ENVIRONMENT="$SOFTWARE_PATH/."

julia --project=$ENVIRONMENT $SOFTWARE
