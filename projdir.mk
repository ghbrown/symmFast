# -*- mode: makefile-gmake -*-
PROJ_DIR_TMP := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROJ_DIR     ?= $(PROJ_DIR_TMP)
