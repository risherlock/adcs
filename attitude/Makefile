CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -Werror -pedantic -lm

TARGET = test
BUILD_DIR = build

INC = \

SRC = \
dcm.c \
maps.c \
quat.c \
test_maps.c

ifeq ($(OS),Windows_NT)
	TARGET := $(TARGET).exe
	M = del /q
	KDIR = mkdir
else
	RM = rm -rf
	MKDIR = mkdir -p
endif

all: $(TARGET)

$(TARGET): $(SRC)
	@if not exist $(BUILD_DIR) $(MKDIR) $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) -o $(BUILD_DIR)/$(TARGET) $(SRC)

run: all
	$(BUILD_DIR)/$(TARGET)

clean:
	$(RM) $(BUILD_DIR)

.PHONY: all clean
