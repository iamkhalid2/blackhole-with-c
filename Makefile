# Black Hole Simulation - Makefile

CC = gcc
CFLAGS = -O2 -Wall -Wextra
LDFLAGS = -lm

SRC_DIR = src
OUT_DIR = output

TARGET = blackhole

SRCS = $(SRC_DIR)/main.c
HEADERS = $(SRC_DIR)/vec3.h $(SRC_DIR)/ray.h $(SRC_DIR)/camera.h \
          $(SRC_DIR)/framebuffer.h $(SRC_DIR)/blackhole.h $(SRC_DIR)/scene.h

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OUT_DIR)/*.ppm
