# Compiler
CC = g++

# Compiler flags
CFLAGS = -std=c++11

# Executable name
EXEC = RayTracer.exe

# Directory containing test files
TEST_DIR = Assignment3-Tests-and-Keys
TEST_FILES = $(wildcard $(TEST_DIR)/*.txt)

# Rule to compile the code
$(EXEC): raytrace.cpp
	$(CC) $(CFLAGS) raytrace.cpp -o $(EXEC)

# Target to run tests
test: $(EXEC)
	@for test_file in $(TEST_FILES); do \
		echo "Running test: $$test_file"; \
		./$(EXEC) $$test_file; \
	done

# Clean generated files
clean:
	rm -f $(EXEC)
