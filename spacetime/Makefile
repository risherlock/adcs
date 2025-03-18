CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -Werror
LDFLAGS = -lm

TARGET = tests
BUILD_DIR = build

INC = \
-I../attitude \
-Iiau06

SRC = \
../attitude/dcm.c \
../attitude/maps.c \
iau06/iau06.c \
iau06/iau06_data.c \
time.c \
frame.c \
test_frame.c

ifeq ($(OS),Windows_NT)
    TARGET := $(TARGET).exe
    RM = del /q
    MKDIR = mkdir
else
    RM = rm -rf
    MKDIR = mkdir
endif

all: $(TARGET)

$(TARGET): $(SRC)
	$(MKDIR) $(BUILD_DIR) || true
	$(CC) $(CFLAGS) $(INC) -o $(BUILD_DIR)/$(TARGET) $(SRC) $(LDFLAGS)

run: all
	$(BUILD_DIR)/$(TARGET)

clean:
	rm -r $(BUILD_DIR) || true

.PHONY: all run clean
