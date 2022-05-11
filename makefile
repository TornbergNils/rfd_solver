SRC_DIR := src
OBJ_DIR := obj
DATA_DIR := data
PLOT_DIR := figures
BIN_DIR := .

CC = g++
CFLAGS := -Wall -g -pedantic -O3
LDLIBS := -lm -g
LDFLAGS := -I./*

EXE := $(BIN_DIR)/program.bin
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

.PHONY: all clean

all: $(EXE) $(DATA_DIR) $(PLOT_DIR)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@
	touch ./src/main.cpp


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR) $(DATA_DIR) $(PLOT_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(OBJ_DIR) $(DATA_DIR) $(EXE) config.csv

-include $(OBJ:.o=.d)
