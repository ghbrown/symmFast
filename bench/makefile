# -*- mode: makefile-gmake -*-
.DELETE_ON_ERROR:

include ../projdir.mk

SRC_EXAMPLES_CXX = $(shell find . -name "*.cpp")

include $(PROJ_DIR)/lib/symmFast/conf/projvariables
include $(PROJ_DIR)/lib/symmFast/conf/projrules
include $(PROJ_DIR)/lib/symmFast/conf/projconf

CLEANABLE = $(patsubst %.cpp,%.dSYM,$(SRC_EXAMPLES_CXX))
CLEANABLE += $(patsubst %.cpp,%,$(SRC_EXAMPLES_CXX))

clean:
	-@rm -rf $(CLEANABLE)
