CXX = g++
CXXFLAGS = -O2 -std=c++11

BIN_DIR = bin
SRC_DIR = src/backend

TARGETS = $(BIN_DIR)/MafftGapConverter $(BIN_DIR)/concatenate

all: $(TARGETS)

$(BIN_DIR)/MafftGapConverter: $(SRC_DIR)/MafftGapConverter.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN_DIR)/concatenate: $(SRC_DIR)/concatenate.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)
