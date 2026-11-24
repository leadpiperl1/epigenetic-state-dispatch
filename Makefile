CC ?= gcc
NM ?= nm
PYTHON ?= python3
JULIA ?= julia
BUILD_DIR ?= $(HOME)/.cache/epigenetic_dispatch_v1_0_2
CFLAGS := -std=c99 -O3 -Wall -Wextra -Wconversion -Werror

.PHONY: all verify clean

all: $(BUILD_DIR)/active_dna_hopf_gate_v1_0_1

$(BUILD_DIR)/active_dna_hopf_gate_v1_0_1: src/active_dna_hopf_gate_v1_0_1.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $< -o $@
	@if $(NM) -u $@ | grep -Eq '\b(malloc|calloc|realloc|free)\b'; then \
		echo "[FAIL] Dynamic heap symbol found in $@" >&2; \
		exit 1; \
	fi
	@echo "[OK] C99 active_dna_hopf_gate_v1_0_1 compiled (malloc = 0)."

verify: all
	$(BUILD_DIR)/active_dna_hopf_gate_v1_0_1
	$(PYTHON) benchmarks/verify_active_dna_mori_zwanzig_v1_0_1.py --seed 42 > $(BUILD_DIR)/active_dna_receipt.json
	$(PYTHON) benchmarks/verify_epigenetic_v1_0_1.py --seed 42
	$(JULIA) --project=julia julia/ActiveDNAWaveguide_v1_0_0.jl
	env -u BASH_ENV bash "$(CURDIR)/scripts/audit_wfx_loot.sh" "$(CURDIR)"
	@echo "[ === EPIGENETIC STATE DISPATCH v1.0.2 VERIFIED (malloc = 0, Julia stdlib, --seed 42, CODE_SHOKUNIN_UNLATCHED = 255) === ]"

clean:
	rm -rf "$(BUILD_DIR)"
