CC = g++ -lm -D__STDC_VERSION__
CFLAGS = -Wall -Wextra -pedantic -fsanitize=undefined,address -lm

TARGET = spg4
BUILD_DIR = build

INC = \
-I../spacetime 

SRC = \
../spacetime/time.c \
SGP4.c \
TLE.cpp \
main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) -o $(BUILD_DIR)/$(TARGET) $(SRC)

run: clean all
	./$(BUILD_DIR)/$(TARGET)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean