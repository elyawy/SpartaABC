CC = g++
 
# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g

# The build target 
TARGET = SpartaABC

all: $(TARGET)

$(TARGET): src/$(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) src/*.cpp

clean:
	$(RM) $(TARGET)
