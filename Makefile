# MetaCSST Makefile
# Bioinformatics tool for DGR prediction using GHMM
# Supports both original C-style and modern C++20 versions

# =============================================================================
# Configuration
# =============================================================================

# Compiler settings
CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -pthread
LDFLAGS := -pthread

# Alternative: Debug build
# CXXFLAGS := -std=c++20 -g -O0 -Wall -Wextra -pthread -DDEBUG

# Source directories
SRCDIR := src
EXAMPLEDIR := example

# Build targets
TARGET_MAIN := MetaCSSTmain
TARGET_SUB := MetaCSSTsub

# Modern C++20 sources
MODERN_MAIN_SRC := $(SRCDIR)/main_modern.cpp
MODERN_SUB_SRC := $(SRCDIR)/sub_modern.cpp
MODERN_HEADERS := $(SRCDIR)/ghmm_modern.h $(SRCDIR)/fun_modern.h

# Original C sources (for reference/comparison)
ORIG_MAIN_SRC := $(SRCDIR)/main.cpp
ORIG_SUB_SRC := $(SRCDIR)/sub.cpp
ORIG_HEADERS := $(SRCDIR)/ghmm.h $(SRCDIR)/fun.h

# Default target
.PHONY: all clean test install modern original help

# =============================================================================
# Default Build (Modern C++20)
# =============================================================================

all: modern

modern: $(TARGET_MAIN) $(TARGET_SUB)

$(TARGET_MAIN): $(MODERN_MAIN_SRC) $(MODERN_HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built: $@ (C++20)"

$(TARGET_SUB): $(MODERN_SUB_SRC) $(MODERN_HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built: $@ (C++20)"

# =============================================================================
# Original Version (for comparison)
# =============================================================================

original: MetaCSSTmain_orig MetaCSSTsub_orig

MetaCSSTmain_orig: $(ORIG_MAIN_SRC) $(ORIG_HEADERS)
	$(CXX) -O2 -lpthread -o $@ $<
	@echo "Built: $@ (original)"

MetaCSSTsub_orig: $(ORIG_SUB_SRC) $(ORIG_HEADERS)
	$(CXX) -O2 -lpthread -o $@ $<
	@echo "Built: $@ (original)"

# =============================================================================
# Testing
# =============================================================================

TESTDIR := test_output
TESTCONFIG := arg.config
TESTINPUT := $(EXAMPLEDIR)/hv29.fa

# Quick test with 1 thread
test: modern
	@echo "Running quick test..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/test_out -thread 1
	@echo "Test output in $(TESTDIR)/test_out/"

# Full test with example data
test-full: modern
	@echo "Running full test pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/raw -thread 4
	perl $(SRCDIR)/callVR.pl $(TESTDIR)/raw/raw.gtf $(TESTINPUT) $(TESTDIR)/tmp1
	perl $(SRCDIR)/removeRepeat.pl $(TESTDIR)/tmp1 $(TESTDIR)/final.gtf
	@echo "Full test complete. Results in $(TESTDIR)/final.gtf"

# Verify output against example
verify: test-full
	@echo "Verifying output..."
	@python3 -c "
	import sys
	with open('$(TESTDIR)/tmp1', 'r') as f:
	    test_lines = set(f.readlines())
	with open('$(EXAMPLEDIR)/out_tmp2.txt', 'r') as f:
	    expected_lines = set(f.readlines())
	common = len(test_lines & expected_lines)
	total = len(expected_lines)
	print(f'Match: {common}/{total} lines ({100*common/total:.1f}%)')
	sys.exit(0 if common == total else 1)
	"

# =============================================================================
# Installation
# =============================================================================

PREFIX ?= /usr/local
BINDIR := $(PREFIX)/bin
DATADIR := $(PREFIX)/share/metacsst

install: modern
	@echo "Installing to $(PREFIX)..."
	@install -d $(BINDIR)
	@install -d $(DATADIR)
	@install -m 755 $(TARGET_MAIN) $(BINDIR)/
	@install -m 755 $(TARGET_SUB) $(BINDIR)/
	@install -m 644 *.config $(DATADIR)/
	@install -m 755 $(SRCDIR)/*.pl $(BINDIR)/
	@echo "Installation complete."
	@echo "Binaries installed to: $(BINDIR)"
	@echo "Data files installed to: $(DATADIR)"

uninstall:
	@echo "Removing from $(PREFIX)..."
	@rm -f $(BINDIR)/$(TARGET_MAIN)
	@rm -f $(BINDIR)/$(TARGET_SUB)
	@rm -f $(BINDIR)/callVR.pl
	@rm -f $(BINDIR)/removeRepeat.pl
	@rm -rf $(DATADIR)
	@echo "Uninstallation complete."

# =============================================================================
# Cleanup
# =============================================================================

clean:
	@echo "Cleaning build artifacts..."
	@rm -f $(TARGET_MAIN) $(TARGET_SUB)
	@rm -f MetaCSSTmain_orig MetaCSSTsub_orig
	@rm -rf $(TESTDIR)
	@rm -f *.o *.a
	@echo "Clean complete."

distclean: clean
	@rm -rf test_baseline test_modern test_fix test_v2 test_sub
	@echo "Deep clean complete."

# =============================================================================
# Development Helpers
# =============================================================================

format:
	@echo "Formatting code..."
	@clang-format -i $(SRCDIR)/*.cpp $(SRCDIR)/*.h 2>/dev/null || \
		echo "clang-format not installed. Skipping format."

check:
	@echo "Running static analysis..."
	@cppcheck --enable=all --suppress=missingIncludeSystem $(SRCDIR) 2>/dev/null || \
		echo "cppcheck not installed. Skipping analysis."

deps:
	@$(CXX) -MM $(MODERN_MAIN_SRC) $(MODERN_SUB_SRC) -I$(SRCDIR) 2>/dev/null || true

# =============================================================================
# Documentation
# =============================================================================

help:
	@echo "MetaCSST Build System"
	@echo "====================="
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build modern C++20 version (default)"
	@echo "  modern       - Build C++20 version with memory safety"
	@echo "  original     - Build original C-style version"
	@echo "  test         - Quick test with 1 thread"
	@echo "  test-full    - Full test with 4 threads and post-processing"
	@echo "  verify       - Compare output against example data"
	@echo "  install      - Install binaries to system (PREFIX=$(PREFIX))"
	@echo "  uninstall    - Remove installed files"
	@echo "  clean        - Remove build artifacts"
	@echo "  distclean    - Deep clean including test outputs"
	@echo "  format       - Format code with clang-format"
	@echo "  check        - Run static analysis"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make                    # Build modern version"
	@echo "  make test-full          # Run complete test"
	@echo "  make install            # Install system-wide"
	@echo "  make clean              # Clean up"

# =============================================================================
# Automatic Dependencies
# =============================================================================

# Include auto-generated dependencies if they exist
-include .deps

.deps: $(MODERN_MAIN_SRC) $(MODERN_SUB_SRC) $(MODERN_HEADERS)
	@$(CXX) -MM $(MODERN_MAIN_SRC) $(MODERN_SUB_SRC) > .deps 2>/dev/null || true

# Mark generated files as secondary (don't delete intermediate files)
.SECONDARY:

# Disable built-in rules for faster builds
.SUFFIXES:
