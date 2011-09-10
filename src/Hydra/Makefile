CXX = g++
CXXFLAGS = -O2
LDFLAGS = -lz 

UTILITIES_DIR = ../common/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES =

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= HydraRouterMain.cpp HydraMain.cpp HydraAssemblerMain.cpp Hydra.cpp Sort.cpp Ancillary.cpp gzstream.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= hydra-ext

all: hydra-ext hydra-assembler hydra-router

.PHONY: all

hydra-router: $(OBJ_DIR)/HydraRouterMain.o $(OBJ_DIR)/Hydra.o $(OBJ_DIR)/Sort.o $(OBJ_DIR)/Ancillary.o $(OBJ_DIR)/gzstream.o
	@echo "  * linking hydra-router"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

hydra-assembler: $(OBJ_DIR)/HydraAssemblerMain.o $(OBJ_DIR)/Hydra.o $(OBJ_DIR)/Sort.o $(OBJ_DIR)/Ancillary.o $(OBJ_DIR)/gzstream.o
	@echo "  * linking hydra-assembler"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

hydra-ext: $(OBJ_DIR)/HydraMain.o $(OBJ_DIR)/Hydra.o $(OBJ_DIR)/Sort.o $(OBJ_DIR)/Ancillary.o $(OBJ_DIR)/gzstream.o
	@echo "  * linking hydra-ext"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)
		
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
