# TelomereHunter Optimized C++ Filter
# Makefile for building the multithreaded telomere read filter

CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -pthread
LDFLAGS = -lhts -lz -lpthread

SRC_DIR = src
BIN_DIR = bin
TARGET = $(BIN_DIR)/telomere_filter

SOURCES = $(SRC_DIR)/telomere_filter.cpp
HEADERS = $(SRC_DIR)/telomere_filter.h

.PHONY: all clean install

all: $(TARGET)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(TARGET): $(SOURCES) $(HEADERS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(LDFLAGS)

clean:
	rm -rf $(BIN_DIR)

install: $(TARGET)
	@echo "Installing telomere_filter to /usr/local/bin (requires sudo)"
	sudo cp $(TARGET) /usr/local/bin/

# Build with debug symbols
debug: CXXFLAGS += -g -DDEBUG
debug: clean $(TARGET)

# Show build information
info:
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Libraries: $(LDFLAGS)"
	@echo "Target: $(TARGET)"
