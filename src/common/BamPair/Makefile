CXX = g++
CXXFLAGS = -O3  -Wall -m64
LDFLAGS = 
OBJ_DIR = ../../../obj/
UTILITIES_DIR = ../

INCLUDES = -I$(UTILITIES_DIR)/BamTools/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= BamPair.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean