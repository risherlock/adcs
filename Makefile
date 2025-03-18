CC = g++ -lm -D__STDC_VERSION__
CFLAGS = -Wall -Wextra -lm

TARGET = main
BUILD_DIR = build

INC = \
-Ispacetime \
-Ispacetime/iau06 \
-Iorbit \
-Iattitude \
-Iigrf

SRC = \
spacetime/time.c \
attitude/dcm.c \
attitude/maps.c \
attitude/quat.c \
spacetime/frame.c \
spacetime/iau06/iau06.c \
spacetime/iau06/iau06_data.c \
orbit/SGP4.c \
orbit/TLE.cpp \
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