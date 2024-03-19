# Compiler
CXX := g++

# Executable name
TARGET := main

# Compiler & Linker flags
INC_DIR := ../include
CXXFLAGS := -I$(INC_DIR) $(shell root-config --cflags)
LDFLAGS := $(shell root-config --glibs)

# Source files
SRCS := $(wildcard ../source/*.cpp) main.cpp

# all: $(TARGET) run

# Compiles and links the program
$(TARGET):
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

# # Run the program
# run:
# 	./$(TARGET)

# Cleans the built objects
clean:
	rm -f $(TARGET)

# Declare phony targets
.PHONY: all run clean