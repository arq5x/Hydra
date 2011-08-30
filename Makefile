# ==========================
# Hydra Makefile
# (c) 2009 Aaron Quinlan
# ==========================

# define our object and binary directories
export OBJ_DIR = obj
export BIN_DIR = bin
export SRC_DIR = src

# define some default flags
export CFLAGS ?= -Wall -O3
export CXXFLAGS ?= $(CFLAGS)
export CXX ?= g++


# define our source subdirectories
SUBDIRS = $(SRC_DIR)/Hydra

all:

	@echo "Building Hydra suite:"
	@echo "========================================================="
	
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done



.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean

