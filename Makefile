# -*- mode: makefile-gmake -*-
.SECONDEXPANSION:
.DELETE_ON_ERROR:

ALL: all

THIS_DIR := $(lastword $(MAKEFILE_LIST))
PROJ_DIR_REL ?= $(patsubst %/,%,$(dir $(THIS_DIR)))
PROJ_DIR ?= $(realpath $(PROJ_DIR_REL))

include $(PROJ_DIR)/lib/symmFast/conf/projvariables
-include $(PROJ_DIR)/lib/symmFast/conf/projrules

CXX = clang++
CXXLINKER = clang++
PROJ_CXXFLAGS += -Wall -fPIC -fvisibility=hidden -fdiagnostics-show-template-tree -fsanitize=address -g -std=gnu++17
PROJ_CXX_DEPFLAGS = -MMD -MP -MF $(PROJ_OBJ_DIR)/$*.d
PROJ_LINKER_FLAGS = $(PROJ_LINKER_FLAGS_BASE) -Wl,-bind_at_load -Wl,-multiply_defined,suppress -Wl,-commons,use_dylibs -Wl,-search_paths_first -Wl,-no_compact_unwind	-fsanitize=address
LDLIBS += -Wl,-rpath,$(PROJ_LIB_DIR) -L$(PROJ_LIB_DIR) -lstdc++

SRC_STRUCTURE := $(shell find $(PROJ_SRC_DIR) -type d)
SRC_TREE := $(addsuffix /*,$(SRC_STRUCTURE))
SRC_LEAVES := $(wildcard $(SRC_TREE))
SRC_ALL = $(filter %.cpp,$(SRC_LEAVES))

INCLUDE_TREE :=	$(shell find $(PROJ_INCLUDE_DIR) -name "*.hpp")
INCLUDES_ALL := $(filter %.hpp,$(INCLUDE_TREE))

OBJ := $(subst $(PROJ_SRC_DIR),$(PROJ_OBJ_DIR),$(SRC_ALL:%.cpp=%.o))
DEPS := $(OBJ:.o=.d)

ifeq ($(V),)
  quiet_HELP := "Use \"$(MAKE) V=1\" to see the verbose compile lines.\n"
  quiet = @printf $(quiet_HELP)$(eval quiet_HELP:=)"  %10s %s\n" "$1$2" "$@"; $($1)
else ifeq ($(V),0)		# Same, but do not print any help
  quiet = @printf "  %10s %s\n" "$1$2" "$@"; $($1)
else				# Show the full command line
  quiet = $($1)
endif

PROJ_COMPILE_CXX = $(call quiet,CXX) $(PROJ_CXXPPFLAGS) $(CXXPPFLAGS) $(PROJ_CXXFLAGS) $(CXXFLAGS) $(PROJ_CXX_DEPFLAGS) -c
PROJ_LINK_CXX = $(call quiet,CXXLINKER) $(PROJ_CXXFLAGS) $(CXXFLAGS) $(PROJ_LINKER_FLAGS) $(LDFLAGS)

LIBS = $(PROJ_LIB_DIR)/$(PROJ_LONG_LIB_NAME)

.PHONY: all info clean all-clean libs all-local TAGS

help:
	@printf "Usage: make [MAKE_OPTIONS] [target] (see 'make --help' for MAKE_OPTIONS)\n"
	@printf ""
	@awk '								\
	{								\
	  if ($$0 ~ /^.PHONY: [a-zA-Z\-\0-9]+$$/) {			\
	    helpCommand = substr($$0, index($$0, ":") + 2);		\
	    if (helpMessage) {						\
	      printf "\033[36m%-20s\033[0m %s\n", helpCommand, helpMessage; \
	      helpMessage = "";						\
	    }								\
	  } else if ($$0 ~ /^[a-zA-Z\-\0-9.]+:/) {			\
	    helpCommand = substr($$0, 0, index($$0, ":"));		\
	    if (helpMessage) {						\
	      printf "\033[36m%-20s\033[0m %s\n", helpCommand, helpMessage; \
	      helpMessage = "";						\
	    }								\
	  } else if ($$0 ~ /^##/) {					\
	    if (helpMessage) {						\
	      helpMessage = helpMessage"\n                     "substr($$0, 3); \
	    } else {							\
	      helpMessage = substr($$0, 3);				\
	    }								\
	  } else {							\
	    if (helpMessage) {						\
	      print "\n                     "helpMessage"\n";		\
	    }								\
	    helpMessage = "";						\
	  }								\
	}'								\
	$(MAKEFILE_LIST)

## -- commonly used --

## delete all build-related files (commonly called dist-clean)
all-clean: clean
	@$(RM) -rf $(PROJ_BUILD_DIR)
	@$(RM) -f $(PROJ_DIR)/make.log

## delete all build files from the current arch directory
clean:
	@$(RM) -rf $(PROJ_OBJ_DIR) $(PROJ_LIB_DIR)

## build the library
all:
	@mkdir -p $(PROJ_ARCH)
	@>$(PROJ_ARCH)/make.log
	@ln -sf $(PROJ_ARCH)/make.log make.log
	+@$(OMAKE_SELF_PRINTDIR) PROJ_DIR=$(PROJ_DIR) PROJ_ARCH=$(PROJ_ARCH) all-local 2>&1 | tee -a ${PROJ_ARCH}/make.log

info:
	-@echo "=========================================="
	-@echo Starting make run on `hostname` at `date +'%a, %d %b %Y %H:%M:%S %z'`
	-@echo Machine characteristics: `uname -a`
	-@echo "-----------------------------------------"
	-@echo "Using PROJ directory: ${PROJ_DIR}"
	-@echo "Using PROJ arch:      ${PROJ_ARCH}"
	-@echo "Using PROJ version:   $(PROJ_VERSION)"
	-@echo "------------------------------------------"
	-@echo "Using CXX: $(shell which $(CXX))"
	-@echo "Using CXXFLAGS: $(PROJ_CXXFLAGS) $(CXXFLAGS)"
	-@echo "Using CXX Version: `$(CXX) --version`"
	-@echo "------------------------------------------"
	-@echo "Using MAKE: $(MAKE)"
	-@echo "Using MAKEFLAGS: -j$(MAKE_NP) -l$(MAKE_LOAD) $(MAKEFLAGS)"
	-@echo "=========================================="

all-local: info $(LIBS)

$(LIBS): $(OBJ) | $$(@D)/.DIR
	$(PROJ_LINK_CXX) -o $@ $^ $(LDLIBS)
	@ln -sf $@ $(PROJ_LIB_DIR)/$(PROJ_LIB_NAME)

$(PROJ_OBJ_DIR)/%.o: $(PROJ_SRC_DIR)/%.cpp | $$(@D)/.DIR
	$(PROJ_COMPILE_CXX) -I$(PROJ_INCLUDE_DIR) $< -o $@

.PRECIOUS: %/.DIR

%/.DIR:
	-@mkdir -p $(@D)
	-@touch $@

TAGS: $(SRC_ALL)
	@$(RM) -f ./TAGS
	@etags $(SRC_ALL)

-include $(DEPS)
