
FC=gfortran
FPP=fypp
FPP_EXTENSION=fypp
F_EXTENSION=f90
#FF= -O3

#source directory variables
SRC_DIR=src
HELLO_DIR=$(SRC_DIR)/hello
SERIAL_DIR=$(SRC_DIR)/serial
PARALLEL_DIR=$(SRC_DIR)/parallel

#build directory variables
BUILD_DIR=build
MOD_DIR=$(BUILD_DIR)/mod
OBJ_DIR=$(BUILD_DIR)/obj
EXE_DIR=$(BUILD_DIR)/exe

#groups of programs
HELLO_PROGS="hello"


all:
	#@$(MAKE) -s serial
	#@$(MAKE) -s parallel

hello:
	$(foreach FPP_FILE, $(HELLO_PROGS), \
		$(FPP) $(HELLO_DIR)/$(FPP_FILE).$(FPP_EXTENSION) $(HELLO_DIR)/$(FPP_FILE).f90;)
	$(foreach F_FILE, $(HELLO_PROGS), \
		$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -o $(OBJ_DIR)/$(F_FILE).o -c $(HELLO_DIR)/$(F_FILE).$(F_EXTENSION))

serial:

parallel:

clean:
	@rm -rf $(BUILD_DIR)
	@mkdir $(BUILD_DIR) $(MOD_DIR) $(OBJ_DIR) $(EXE_DIR)
	@find ./ -name '*.$(F_EXTENSION)' -delete
