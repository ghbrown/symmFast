
FC=gfortran
FPP=fypp
FPP_EXTENSION=fypp
F_EXTENSION=f90
#FF= -O3

#source directory variables
SRC_DIR=src
HELLO_DIR=$(SRC_DIR)/hello
CANONICAL_DIR=$(SRC_DIR)/canonical
FAST_DIR=$(SRC_DIR)/fast

#build directory variables
BUILD_DIR=build
MOD_DIR=$(BUILD_DIR)/mod
OBJ_DIR=$(BUILD_DIR)/obj
EXE_DIR=$(BUILD_DIR)/exe

#groups of programs
HELLO_PROGS="hello"
CANONICAL_PROGS="c_d_s_symv"


all:
	#@$(MAKE) -s serial
	#@$(MAKE) -s parallel

kinds:
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -o $(OBJ_DIR)/stdlib_kinds.o -c $(SRC_DIR)/stdlib_kinds.f90

hello:
	$(foreach FPP_FILE, $(HELLO_PROGS), \
		$(FPP) $(HELLO_DIR)/$(FPP_FILE).$(FPP_EXTENSION) $(HELLO_DIR)/$(FPP_FILE).f90;)
	$(foreach F_FILE, $(HELLO_PROGS), \
		$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -o $(OBJ_DIR)/$(F_FILE).o -c $(HELLO_DIR)/$(F_FILE).$(F_EXTENSION))

canonical:
	@$(MAKE) -s kinds
	$(foreach FPP_FILE, $(CANONICAL_PROGS), \
		$(FPP) $(CANONICAL_DIR)/$(FPP_FILE).$(FPP_EXTENSION) $(CANONICAL_DIR)/$(FPP_FILE).f90;)
	$(foreach F_FILE, $(CANONICAL_PROGS), \
		$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -o $(OBJ_DIR)/$(F_FILE).o -c $(CANONICAL_DIR)/$(F_FILE).$(F_EXTENSION))

fast:

clean:
	@rm -rf $(BUILD_DIR)
	@mkdir $(BUILD_DIR) $(MOD_DIR) $(OBJ_DIR) $(EXE_DIR)
	@find $(CANONICAL_DIR) -name '*.$(F_EXTENSION)' -delete
	@find $(FAST_DIR) -name '*.$(F_EXTENSION)' -delete
