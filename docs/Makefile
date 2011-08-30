CXX = g++
CXXFLAGS = -O3
LDFLAGS =

SRC_DIR = ../../src/
UTILITIES_DIR = ../common/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(SRC_DIR)/ExternalMergeSort/ -I$(UTILITIES_DIR)/Utilities/  -I$(UTILITIES_DIR)/Hydra/ -I$(UTILITIES_DIR)/Utilities/gzstream/---

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= HydraMain.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS= Hydra.o Sort.o Ancillary.o gzstream.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= hydra

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)
	
$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/Hydra/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/Utilities/gzstream/
			
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
