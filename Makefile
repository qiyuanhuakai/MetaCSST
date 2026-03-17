# MetaCSST Makefile (modern-only)

CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -pthread
LDFLAGS := -pthread

SRCDIR := src
EXAMPLEDIR := example

TARGET_MAIN := MetaCSSTmain
TARGET_SUB := MetaCSSTsub

MAIN_SRC := $(SRCDIR)/main_modern.cpp
SUB_SRC := $(SRCDIR)/sub_modern.cpp
HEADERS := $(SRCDIR)/ghmm_modern.hpp $(SRCDIR)/fun_modern.hpp $(SRCDIR)/config_modern.hpp

TESTDIR := test_output
TESTCONFIG := config.json
TESTINPUT := $(EXAMPLEDIR)/hv29.fa

PREFIX ?= /usr/local
BINDIR := $(PREFIX)/bin
DATADIR := $(PREFIX)/share/metacsst

.PHONY: all modern test test-full verify verify-toml verify-yaml install uninstall clean help

all: modern

modern: $(TARGET_MAIN) $(TARGET_SUB)

$(TARGET_MAIN): $(MAIN_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built: $@"

$(TARGET_SUB): $(SUB_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built: $@"

test: modern
	@echo "Running quick test..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/quick_raw -thread 1
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/quick_raw/raw.gtf $(TESTINPUT) $(TESTDIR)/quick_final.gtf
	@echo "Quick test output: $(TESTDIR)/quick_final.gtf"

test-full: modern
	@echo "Running full test pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/raw -thread 4
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/raw/raw.gtf $(TESTINPUT) $(TESTDIR)/final.gtf
	@echo "Full test complete. Results in $(TESTDIR)/final.gtf"

verify: test-full
	@echo "Verifying output..."
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'Match: {common}/{total} lines ({100*common/total:.1f}%)'); sys.exit(0 if test_lines==expected_lines else 1)"

verify-toml: modern
	@echo "Verifying TOML pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build config.toml -in $(TESTINPUT) -out $(TESTDIR)/toml_raw -thread 1
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/toml_raw/raw.gtf $(TESTINPUT) $(TESTDIR)/toml_final.gtf
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/toml_final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'TOML Match: {common}/{total} lines ({100*common/total:.1f}%)'); sys.exit(0 if test_lines==expected_lines else 1)"

verify-yaml: modern
	@echo "Verifying YAML pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build config.yaml -in $(TESTINPUT) -out $(TESTDIR)/yaml_raw -thread 1
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/yaml_raw/raw.gtf $(TESTINPUT) $(TESTDIR)/yaml_final.gtf
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/yaml_final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'YAML Match: {common}/{total} lines ({100*common/total:.1f}%)'); sys.exit(0 if test_lines==expected_lines else 1)"

install: modern
	@echo "Installing to $(PREFIX)..."
	@install -d $(BINDIR)
	@install -d $(DATADIR)
	@install -m 755 $(TARGET_MAIN) $(BINDIR)/
	@install -m 755 $(TARGET_SUB) $(BINDIR)/
	@install -m 755 $(SRCDIR)/call_vr.py $(BINDIR)/
	@install -m 644 config.json config.toml config.yaml $(DATADIR)/
	@echo "Installation complete."

uninstall:
	@echo "Removing from $(PREFIX)..."
	@rm -f $(BINDIR)/$(TARGET_MAIN) $(BINDIR)/$(TARGET_SUB) $(BINDIR)/call_vr.py
	@rm -rf $(DATADIR)
	@echo "Uninstallation complete."

clean:
	@echo "Cleaning build artifacts..."
	@rm -f $(TARGET_MAIN) $(TARGET_SUB) .deps
	@rm -rf $(TESTDIR)
	@echo "Clean complete."

help:
	@echo "Targets: all modern test test-full verify verify-toml verify-yaml install uninstall clean"
