# Compiler and flags
CXX = g++
CXXFLAGS = -O3 -std=c++11 -Iinclude
SRC_DIR = src
INCLUDE_DIR = include
OUTPUT_DIR = outputs

# All source files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Executable
TARGET = rasterizer

# Build the project
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile each .cpp file into .o
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean the build
clean:
	rm -f $(TARGET) $(SRC_DIR)/*.o

# Move PPM files to outputs directory
move_ppm:
	mv *.ppm $(OUTPUT_DIR)/

# Run the rasterizer and move PPM files
run: $(TARGET)
	@if [ -z "$(input)" ]; then \
		echo "Usage: make run input=<input_file>"; \
		exit 1; \
	fi; \
	./$(TARGET) $(input); \
	make move_ppm

# Test all XML files in the inputs directory
test: $(TARGET)
	@for file in $(wildcard inputs/*.xml); do \
		echo "Running on $$file..."; \
		make run input=$$file; \
	done