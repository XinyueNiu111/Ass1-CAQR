CC     = gcc
CFLAGS = -O2 -Wall
LIBS   = -framework Accelerate

TARGET = tsqr
SRCS   = tsqr.c

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)