SRC_DIR := src
OBJ_DIR := obj
DATA_DIR := data data/time
PLOT_DIR := figures
BIN_DIR := .

CC = g++
CFLAGS := -Wall -g -pedantic
LDLIBS := -lm
LDFLAGS := -I./*

EXE := $(BIN_DIR)/program.bin
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

.PHONY: all clean

all: clean $(EXE) $(DATA_DIR) $(PLOT_DIR)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR) $(DATA_DIR) $(PLOT_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(OBJ_DIR) $(DATA_DIR) $(EXE)

-include $(OBJ:.o=.d)
