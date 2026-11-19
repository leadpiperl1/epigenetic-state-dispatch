CC ?= gcc
PYTHON ?= python3
BUILD_DIR ?= $(HOME)/.cache/epigenetic_dispatch_v1_0_1
CFLAGS := -std=c99 -O3 -Wall -Wextra -Wconversion -Werror

.PHONY: all verify clean

all: $(BUILD_DIR)/active_dna_hopf_gate_v1_0_1

$(BUILD_DIR)/active_dna_hopf_gate_v1_0_1: src/active_dna_hopf_gate_v1_0_1.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $< -o $@

verify: all
	$(BUILD_DIR)/active_dna_hopf_gate_v1_0_1
	$(PYTHON) benchmarks/verify_active_dna_mori_zwanzig_v1_0_1.py --seed 42 > $(BUILD_DIR)/active_dna_receipt.json
	$(PYTHON) benchmarks/verify_epigenetic_v1_0_1.py --seed 42
	env -u BASH_ENV bash "$(CURDIR)/scripts/audit_wfx_loot.sh" "$(CURDIR)"

clean:
	rm -rf "$(BUILD_DIR)"
