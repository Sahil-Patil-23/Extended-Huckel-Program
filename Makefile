# GNU C++ Compiler
CPP  =  g++

CPP_FLAGS  =  -std=c++14 -g -o

SRC_DIR  =  src
BIN_DIR  = bin

all:
	cd $(SRC_DIR) && $(CPP) $(CPP_FLAGS) ../$(BIN_DIR)/main main.cpp -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo
	cd $(BIN_DIR) && ./main