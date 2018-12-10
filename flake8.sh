#!/bin/bash

clear
flake8 --hang-closing --ignore=W293 $1
