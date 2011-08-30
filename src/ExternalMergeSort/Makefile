OBJ_DIR = ../../../obj
BIN_DIR = ../../../bin

INCLUDES = 

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES = HydraMain.cpp \
          Hydra.cpp \
          Ancilary.cpp \
          Sort.cpp \
          SequenceUtilities.cpp \
          gzstream.cpp \

OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
EXT_OBJECTS = $(OBJ_DIR)/*.o
PROGRAM = hydra

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
@echo " * linking $(PROGRAM)"
@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
@echo " * compiling" $(*F).cpp
@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)


clean:
@echo "Cleaning up."
@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean