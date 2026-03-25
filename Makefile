# MetaCSST Makefile (modern-only)

CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -pthread
LDFLAGS := -pthread -lz

SRCDIR := src
EXAMPLEDIR := example

# Third-party dependency management (deterministic local vendor)
THIRD_PARTY_DIR := $(SRCDIR)/third_party
VENDOR_DIR := $(THIRD_PARTY_DIR)/vendor
BUILD_DIR := $(THIRD_PARTY_DIR)/build
STAMP_DIR := $(THIRD_PARTY_DIR)/stamps

NLOHMANN_JSON_VERSION := v3.11.3
NLOHMANN_JSON_ARCHIVE := $(VENDOR_DIR)/nlohmann-json-$(NLOHMANN_JSON_VERSION).tar.xz
NLOHMANN_JSON_DIR := $(VENDOR_DIR)/json
NLOHMANN_JSON_INC := $(NLOHMANN_JSON_DIR)/single_include

TOMLPP_VERSION := v3.4.0
TOMLPP_ARCHIVE := $(VENDOR_DIR)/tomlplusplus-$(TOMLPP_VERSION).tar.gz
TOMLPP_DIR := $(VENDOR_DIR)/tomlplusplus
TOMLPP_INC := $(TOMLPP_DIR)/include

YAML_CPP_VERSION := 0.8.0
YAML_CPP_TAG := yaml-cpp-$(YAML_CPP_VERSION)
YAML_CPP_ARCHIVE := $(VENDOR_DIR)/yaml-cpp-$(YAML_CPP_VERSION).tar.gz
YAML_CPP_SRC_DIR := $(VENDOR_DIR)/yaml-cpp-$(YAML_CPP_VERSION)
YAML_CPP_BUILD_DIR := $(BUILD_DIR)/yaml-cpp
YAML_CPP_LIB := $(YAML_CPP_BUILD_DIR)/libyaml-cpp.a
YAML_CPP_INC := $(YAML_CPP_SRC_DIR)/include

THIRD_PARTY_CFLAGS := -I$(NLOHMANN_JSON_INC) -I$(TOMLPP_INC) -I$(YAML_CPP_INC)
THIRD_PARTY_LIBS := $(YAML_CPP_LIB)

CXXFLAGS += $(THIRD_PARTY_CFLAGS)
LDFLAGS += $(THIRD_PARTY_LIBS)

TARGET_MAIN := MetaCSSTmain
TARGET_SUB := MetaCSSTsub

MAIN_SRC := $(SRCDIR)/main_modern.cpp
SUB_SRC := $(SRCDIR)/sub_modern.cpp
HEADERS := $(SRCDIR)/ghmm_modern.hpp $(SRCDIR)/fun_modern.hpp $(SRCDIR)/config_modern.hpp \
	$(SRCDIR)/app_common.hpp $(SRCDIR)/fasta_runtime.hpp $(SRCDIR)/thread_runtime.hpp \
	$(SRCDIR)/main_scan.hpp $(SRCDIR)/sub_scan.hpp $(SRCDIR)/scan_pipeline.hpp \
	$(SRCDIR)/main_formatter.hpp $(SRCDIR)/sub_formatter.hpp

TESTDIR := test_output
TESTCONFIG := config.json
TESTINPUT := $(EXAMPLEDIR)/hv29.fa
TESTINPUT_GZ := $(EXAMPLEDIR)/hv29.fa.gz

PREFIX ?= /usr/local
BINDIR := $(PREFIX)/bin
DATADIR := $(PREFIX)/share/metacsst

.PHONY: deps modern test verify verify-json verify-toml verify-yaml verify-compressed verify-thread-consistency verify-sub-consistency install uninstall clean help example


deps: $(STAMP_DIR)/deps.ready

$(STAMP_DIR):
	@mkdir -p $(STAMP_DIR)

$(VENDOR_DIR):
	@mkdir -p $(VENDOR_DIR)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(NLOHMANN_JSON_ARCHIVE): | $(VENDOR_DIR)
	@echo "Fetching nlohmann/json $(NLOHMANN_JSON_VERSION)..."
	@curl -L --fail --retry 3 --retry-delay 2 -o $@ https://github.com/nlohmann/json/releases/download/$(NLOHMANN_JSON_VERSION)/json.tar.xz

$(NLOHMANN_JSON_INC): $(NLOHMANN_JSON_ARCHIVE)
	@echo "Extracting nlohmann/json..."
	@rm -rf $(NLOHMANN_JSON_DIR)
	@tar -xf $(NLOHMANN_JSON_ARCHIVE) -C $(VENDOR_DIR)

$(TOMLPP_ARCHIVE): | $(VENDOR_DIR)
	@echo "Fetching toml++ $(TOMLPP_VERSION)..."
	@curl -L --fail --retry 3 --retry-delay 2 -o $@ https://github.com/marzer/tomlplusplus/archive/refs/tags/$(TOMLPP_VERSION).tar.gz

$(TOMLPP_INC): $(TOMLPP_ARCHIVE)
	@echo "Extracting toml++..."
	@rm -rf $(TOMLPP_DIR)
	@tar -xzf $(TOMLPP_ARCHIVE) -C $(VENDOR_DIR)
	@mv $(VENDOR_DIR)/tomlplusplus-3.4.0 $(TOMLPP_DIR)

$(YAML_CPP_ARCHIVE): | $(VENDOR_DIR)
	@echo "Fetching yaml-cpp $(YAML_CPP_TAG)..."
	@curl -L --fail --retry 3 --retry-delay 2 -o $@ https://github.com/jbeder/yaml-cpp/archive/refs/tags/$(YAML_CPP_VERSION).tar.gz

$(YAML_CPP_SRC_DIR): $(YAML_CPP_ARCHIVE)
	@echo "Extracting yaml-cpp..."
	@rm -rf $(YAML_CPP_SRC_DIR)
	@tar -xzf $(YAML_CPP_ARCHIVE) -C $(VENDOR_DIR)
	@test -d $(YAML_CPP_SRC_DIR) || (echo "Error: yaml-cpp source extraction failed at $(YAML_CPP_SRC_DIR)" && exit 1)

$(YAML_CPP_LIB): $(YAML_CPP_SRC_DIR) | $(BUILD_DIR)
	@echo "Building yaml-cpp static library..."
	@command -v cmake >/dev/null 2>&1 || (echo "Error: cmake not found in PATH." && exit 1)
	@test -d $(YAML_CPP_SRC_DIR) || (echo "Error: yaml-cpp source directory not found: $(YAML_CPP_SRC_DIR)" && exit 1)
	@rm -rf $(YAML_CPP_BUILD_DIR)
	@mkdir -p $(YAML_CPP_BUILD_DIR)
	@cmake -E chdir $(YAML_CPP_BUILD_DIR) cmake $(abspath $(YAML_CPP_SRC_DIR)) -DCMAKE_BUILD_TYPE=Release -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_BUILD_SHARED_LIBS=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5
	@cmake --build $(YAML_CPP_BUILD_DIR) 

$(STAMP_DIR)/deps.ready: $(NLOHMANN_JSON_INC) $(TOMLPP_INC) $(YAML_CPP_LIB) | $(STAMP_DIR)
	@touch $@

modern: deps $(TARGET_MAIN) $(TARGET_SUB)

$(TARGET_MAIN): $(MAIN_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built: $@"

$(TARGET_SUB): $(SUB_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built: $@"


verify: verify-json verify-toml verify-yaml verify-compressed verify-thread-consistency verify-sub-consistency
	@echo "All verify pipelines passed."

.NOTPARALLEL: verify verify-json verify-toml verify-yaml verify-compressed verify-thread-consistency verify-sub-consistency

test: verify
	@echo "test target is merged into verify (alias)."

verify-json: modern
	@echo "Running full test pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/raw -thread 4
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/raw/raw.gtf $(TESTINPUT) $(TESTDIR)/final.gtf
	@echo "Full test complete. Results in $(TESTDIR)/final.gtf"
	@echo "Verifying JSON pipeline..."
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'JSON Match: {common}/{total} lines ({100*common/total:.1f}%) [strict line-content multiset equality, order-insensitive]'); sys.exit(0 if test_lines==expected_lines else 1)"

verify-toml: modern
	@echo "Verifying TOML pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build config.toml -in $(TESTINPUT) -out $(TESTDIR)/toml_raw -thread 1
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/toml_raw/raw.gtf $(TESTINPUT) $(TESTDIR)/toml_final.gtf
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/toml_final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'TOML Match: {common}/{total} lines ({100*common/total:.1f}%) [strict line-content multiset equality, order-insensitive]'); sys.exit(0 if test_lines==expected_lines else 1)"

verify-yaml: modern
	@echo "Verifying YAML pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build config.yaml -in $(TESTINPUT) -out $(TESTDIR)/yaml_raw -thread 1
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/yaml_raw/raw.gtf $(TESTINPUT) $(TESTDIR)/yaml_final.gtf
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/yaml_final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'YAML Match: {common}/{total} lines ({100*common/total:.1f}%) [strict line-content multiset equality, order-insensitive]'); sys.exit(0 if test_lines==expected_lines else 1)"

verify-compressed: modern
	@echo "Running compressed FASTA pipeline..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT_GZ) -out $(TESTDIR)/gz_raw -thread 1
	python3 $(SRCDIR)/call_vr.py $(TESTDIR)/gz_raw/raw.gtf $(TESTINPUT) $(TESTDIR)/gz_final.gtf
	@echo "Compressed pipeline output: $(TESTDIR)/gz_final.gtf"
	@echo "Verifying compressed FASTA pipeline..."
	@python3 -c "import sys; from collections import Counter; test_lines=Counter(open('$(TESTDIR)/gz_final.gtf','r').readlines()); expected_lines=Counter(open('$(EXAMPLEDIR)/out-DGR.gtf','r').readlines()); common=sum((test_lines & expected_lines).values()); total=sum(expected_lines.values()); print(f'GZ Match: {common}/{total} lines ({100*common/total:.1f}%) [strict line-content multiset equality, order-insensitive]'); sys.exit(0 if test_lines==expected_lines else 1)"

verify-thread-consistency: modern
	@echo "Verifying main thread consistency (1 vs 4)..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/raw_t1 -thread 1
	./$(TARGET_MAIN) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/raw_t4 -thread 4
	@python3 -c "import sys; from collections import Counter; a=Counter(open('$(TESTDIR)/raw_t1/raw.gtf','r').readlines()); b=Counter(open('$(TESTDIR)/raw_t4/raw.gtf','r').readlines()); common=sum((a & b).values()); total=max(sum(a.values()), sum(b.values())); print(f'MAIN Thread Consistency: {common}/{total} lines ({100*common/total if total else 100:.1f}%) [strict line-content multiset equality, order-insensitive]'); sys.exit(0 if a==b else 1)"

verify-sub-consistency: modern
	@echo "Verifying sub thread consistency (1 vs 2)..."
	@mkdir -p $(TESTDIR)
	./$(TARGET_SUB) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/sub_t1 -thread 1
	./$(TARGET_SUB) -build $(TESTCONFIG) -in $(TESTINPUT) -out $(TESTDIR)/sub_t2 -thread 2
	@python3 -c "import sys; from collections import Counter; a=Counter(open('$(TESTDIR)/sub_t1/out.txt','r').readlines()); b=Counter(open('$(TESTDIR)/sub_t2/out.txt','r').readlines()); common=sum((a & b).values()); total=max(sum(a.values()), sum(b.values())); print(f'SUB Thread Consistency: {common}/{total} lines ({100*common/total if total else 100:.1f}%) [strict line-content multiset equality, order-insensitive]'); sys.exit(0 if a==b else 1)"




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
	@rm -rf $(THIRD_PARTY_DIR)
	@echo "Clean complete."

example:
	@echo "MetaCSST usage examples:"
	@echo "  1) Build models + run full prediction (JSON config)"
	@echo "     ./MetaCSSTmain -build config.json -in example/hv29.fa -out example/run_out -thread 4"
	@echo "  2) Build models + run with TOML config"
	@echo "     ./MetaCSSTmain -build config.toml -in example/hv29.fa -out example/run_out_toml -thread 2"
	@echo "  3) Build models + run with YAML config"
	@echo "     ./MetaCSSTmain -build config.yaml -in example/hv29.fa -out example/run_out_yaml -thread 2"
	@echo "  4) Sub-structure scan only (TR/VR/RT)"
	@echo "     ./MetaCSSTsub -build config.json -in example/hv29.fa -out example/sub_out -thread 2"
	@echo "  5) Post-process raw output"
	@echo "     python3 src/call_vr.py example/run_out/raw.gtf example/hv29.fa example/final.gtf"

help:
	@echo "MetaCSST make targets:"
	@echo "  deps       : Fetch/build local third-party parser dependencies"
	@echo "  modern     : Build MetaCSSTmain and MetaCSSTsub"
	@echo "  verify     : Build + run all verification pipelines"
	@echo "  example    : Print runnable usage examples"
	@echo "  install    : Install binaries/configs"
	@echo "  uninstall  : Remove installed files"
	@echo "  clean      : Remove binaries and test_output"
	@echo ""
	@echo "Verify criterion: strict line-content multiset equality (order-insensitive)."
