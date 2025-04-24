CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -Werror -lm -O3

TARGET = test
BUILD_DIR = build

INC =
SRC = dcm.c euler.c quat.c axan.c test_att.c

OBJ = $(SRC:%.c=$(BUILD_DIR)/%.o)

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(INC) -o $(BUILD_DIR)/$(TARGET) $(OBJ)
