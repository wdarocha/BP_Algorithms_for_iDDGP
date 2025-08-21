# ====================================================================
# BP_Algorithms_for_iDDGP -- Build System (Linux/macOS/Windows via Git Bash)
#
# Layout:
#   include/adt/*.h		public ADT headers
#   include/algorithms/*.h	public algorithm headers
#   src/adt/**.c		→ build/lib/libadt.a
#   src/algorithms/**.c		→ build/lib/libalgorithms.a (includes common.c + subfolders)
#   src/algorithms/{iBP,iTBP}/	per-algo code + version files
#   src/main.c			→ build/bin/main
#
# Project headers: use #include "adt/..." and #include "algorithms/..."
# Only -Iinclude is passed.
#
# Outputs:
#   build/obj/**.{o,d}	objects + depfiles
#   build/lib/*.a	static libraries
#   build/bin/main      executable
#   build/generated/*.h generated headers (project_version.h)
#
# Versioning:
#   Project: VERSION at repo root
#   iBP:  src/algorithms/iBP/	→ IBP_VERSION
#   iTBP: src/algorithms/iTBP/	→ ITBP_VERSION
#   *_RELEASE_DATE missing	→ use mtime of the corresponding *_VERSION file;
#   otherwise fallback to project date. (Cross-platform stat/date)
# ====================================================================

SHELL := /bin/bash
MAKEFLAGS += --no-print-directory
CC := gcc
AR := ar

# ---------------- OS detection ----------------
UNAME_S := $(shell uname -s 2>/dev/null || echo Unknown)

# Cross-platform: ISO date (YYYY-MM-DD) from file mtime
define mtime_iso
$(shell \
	if [ -f $(1) ]; then \
		if [ "$(UNAME_S)" = "Darwin" ]; then date -r "$$\(stat -f %m $(1)\)" +%F; \
		else date -d "@$$\(stat -c %Y $(1)\)" +%F 2>/dev/null || date -u -r "$$\(stat -c %Y $(1)\)" +%F 2>/dev/null; \
		fi; \
	fi \
)
endef

# ---------------- Dirs ----------------
INC_DIR := include
SRC_DIR := src
ADT_DIR := $(SRC_DIR)/adt
ALG_DIR := $(SRC_DIR)/algorithms
IBP_DIR := $(ALG_DIR)/iBP
ITP_DIR := $(ALG_DIR)/iTBP

MAIN_SRC := $(SRC_DIR)/main.c

BUILD_DIR := build
OBJ_DIR   := $(BUILD_DIR)/obj
LIB_DIR   := $(BUILD_DIR)/lib
BIN_DIR   := $(BUILD_DIR)/bin
GEN_DIR   := $(BUILD_DIR)/generated

.PHONY: dirs
dirs:
	@mkdir -p $(LIB_DIR) $(BIN_DIR) $(GEN_DIR) $(OBJ_DIR)

# ---------------- Project versioning ----------------
VERSION  := $(shell [ -f VERSION ] && cat VERSION || echo 0.0.0)
GIT_HASH := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GIT_DESC := $(shell git describe --tags --dirty --always 2>/dev/null || echo $(GIT_HASH))

RELEASE_DATE := $(shell \
	if [ -f RELEASE_DATE ]; then cat RELEASE_DATE; \
	elif [ -f VERSION ]; then echo "$(call mtime_iso,VERSION)"; \
	elif git rev-parse --git-dir >/dev/null 2>&1; then \
		D=$$(git for-each-ref refs/tags/v$(VERSION) --format='%(taggerdate:short)' | head -n1); \
		if [ -z "$$D" ]; then D=$$(git log -1 --date=short --format=%cd 2>/dev/null); \
		fi; \
		echo $$D; \
	else date +%F; \
	fi \
)

IBP_VERSION  := $(shell [ -f $(IBP_DIR)/IBP_VERSION ] && cat $(IBP_DIR)/IBP_VERSION || echo $(VERSION))
ITBP_VERSION := $(shell [ -f $(ITP_DIR)/ITBP_VERSION ] && cat $(ITP_DIR)/ITBP_VERSION || echo $(VERSION))

IBP_RELEASE_DATE := $(shell \
	if [ -f $(IBP_DIR)/IBP_RELEASE_DATE ]; then cat $(IBP_DIR)/IBP_RELEASE_DATE; \
	elif [ -f $(IBP_DIR)/IBP_VERSION ]; then echo "$(call mtime_iso,$(IBP_DIR)/IBP_VERSION)"; \
	else echo $(RELEASE_DATE); \
	fi \
)

ITBP_RELEASE_DATE := $(shell \
	if [ -f $(ITP_DIR)/ITBP_RELEASE_DATE ]; then cat $(ITP_DIR)/ITBP_RELEASE_DATE; \
	elif [ -f $(ITP_DIR)/ITBP_VERSION ]; then echo "$(call mtime_iso,$(ITP_DIR)/ITBP_VERSION)"; \
	else echo $(RELEASE_DATE); \
	fi \
)

# Generated header
VER_HDR  := $(GEN_DIR)/project_version.h
VER_DEPS := $(wildcard VERSION RELEASE_DATE \
            $(IBP_DIR)/IBP_VERSION $(IBP_DIR)/IBP_RELEASE_DATE \
            $(ITP_DIR)/ITBP_VERSION $(ITP_DIR)/ITBP_RELEASE_DATE)
CPPFLAGS += -I$(GEN_DIR)

$(VER_HDR): Makefile $(VER_DEPS) | dirs
	@{ \
	  echo '#pragma once'; \
	  echo '#define PROJECT_VERSION        "$(VERSION)"'; \
	  echo '#define PROJECT_GIT_HASH       "$(GIT_HASH)"'; \
	  echo '#define PROJECT_GIT_DESCRIBE   "$(GIT_DESC)"'; \
	  echo '#define PROJECT_RELEASE_DATE   "$(RELEASE_DATE)"'; \
	  echo ''; \
	  echo '/* Per-algorithm versioning */'; \
	  echo '#define IBP_VERSION        "$(IBP_VERSION)"'; \
	  echo '#define IBP_RELEASE_DATE   "$(IBP_RELEASE_DATE)"'; \
	  echo '#define ITBP_VERSION       "$(ITBP_VERSION)"'; \
	  echo '#define ITBP_RELEASE_DATE  "$(ITBP_RELEASE_DATE)"'; \
	} > $(VER_HDR)
	@echo "[gen] $(VER_HDR)"

# ---------------- Outputs ----------------
LIBADT := $(LIB_DIR)/libadt.a
LIBALG := $(LIB_DIR)/libalgorithms.a
BIN := $(BIN_DIR)/main
MAPFILE := $(BUILD_DIR)/main.map

# ---------------- Build mode ----------------
MODE ?= release
WARN := -Wall -Wextra -Wpedantic -Wshadow \
        -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls \
        -Werror=implicit-function-declaration

CPPFLAGS += -I$(INC_DIR) -MMD -MP

ifeq ($(MODE), release)
	CFLAGS  := -std=c11 -O3 -march=native -flto -ffunction-sections -fdata-sections $(WARN)
	LDFLAGS := -flto -Wl,--gc-sections -Wl,-O2 -Wl,-Map,$(MAPFILE)
else
	CFLAGS  := -std=c11 -O0 -g -fno-inline -fno-omit-frame-pointer $(WARN)
	LDFLAGS := -Wl,-Map,$(MAPFILE)
endif

# ---------------- Scientific libs ----------------
LAPACKE_LIBS := $(shell pkg-config --silence-errors --libs lapacke)

ifeq ($(strip $(LAPACKE_LIBS)),)
	LDLIBS  := -llapacke -llapack -lblas
	LDFLAGS += -Wl,--no-as-needed
else
	LDLIBS  := $(LAPACKE_LIBS)
endif

MATHLIB := -lm

# ---------------- Sources ----------------
normaliz = sed 's|^\./||'
ADT_SRCS  := $(shell find $(ADT_DIR) -type f -name '*.c' -print | $(normaliz))
ALG_SRCS  := $(shell find $(ALG_DIR) -type f -name '*.c' -print | $(normaliz))
MAIN_SRCS := $(MAIN_SRC)

ADT_OBJS  := $(patsubst $(SRC_DIR)/%, $(OBJ_DIR)/%,$(ADT_SRCS:.c=.o))
ALG_OBJS  := $(patsubst $(SRC_DIR)/%, $(OBJ_DIR)/%,$(ALG_SRCS:.c=.o))
MAIN_OBJS := $(patsubst $(SRC_DIR)/%, $(OBJ_DIR)/%,$(MAIN_SRCS:.c=.o))

DEPS := $(ADT_OBJS:.o=.d) $(ALG_OBJS:.o=.d) $(MAIN_OBJS:.o=.d)

# ---------------- Echo helpers ----------------
ECHO_CC = @echo "[cc] $<"
ECHO_AR = @echo "[ar] $@"
ECHO_LD = @echo "[ld] $@"

# ---------------- Default goals ----------------
.DEFAULT_GOAL := build

.PHONY: build
build: $(BIN) $(LIBALG) $(LIBADT) $(VER_HDR)
	@if [ -n "$$(find $^ -newer .build_stamp 2>/dev/null)" ]; then \
		echo "[build] build updated."; \
		touch .build_stamp; \
	else \
		echo "[build] everything is up to date."; \
	fi

# ---------------- Libraries ----------------
$(LIBADT): $(ADT_OBJS) | dirs
	$(ECHO_AR)
	@$(AR) rcs $@ $^

$(LIBALG): $(ALG_OBJS) | dirs
	$(ECHO_AR)
	@$(AR) rcs $@ $^

# ---------------- Executable ----------------
$(BIN): $(MAIN_OBJS) $(LIBALG) $(LIBADT) | dirs
	$(ECHO_LD)
	@$(CC) $(LDFLAGS) $(MAIN_OBJS) -Wl,--start-group $(LIBALG) $(LIBADT) $(LDLIBS) -Wl,--end-group $(MATHLIB) -o $@
	@echo "[ok] linked $(BIN)"

# ---------------- Compile ----------------
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(VER_HDR) | dirs
	@mkdir -p $(dir $@)
	$(ECHO_CC)
	@$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
	@echo "[ok] compiled $@"

# ---------------- Utilities ----------------
.PHONY: libs clean distclean print-% version help init-algo-meta

libs: $(LIBADT) $(LIBALG)
	@echo "[libs] built: $(LIBADT) and $(LIBALG)"

clean:
	@echo "[clean] removing $(BUILD_DIR) folder and all its subfolders and files"
	@rm -rf $(BUILD_DIR)
	@echo "[clean] done"

print-%:
	@echo "$*=$($*)"

version:
	@echo "Project : $(VERSION) (date: $(RELEASE_DATE))"
	@echo "Git     : $(GIT_HASH) / $(GIT_DESC)"
	@echo "iBP     : $(IBP_VERSION) (date: $(IBP_RELEASE_DATE))"
	@echo "iTBP    : $(ITBP_VERSION) (date: $(ITBP_RELEASE_DATE))"
	@echo "Header  : $(VER_HDR)"

init-algo-meta:
	@mkdir -p $(IBP_DIR) $(ITP_DIR)
	@echo "$(VERSION)" > $(IBP_DIR)/IBP_VERSION
	@echo "$(VERSION)" > $(ITP_DIR)/ITBP_VERSION
	@{ \
		if [ "$(UNAME_S)" = "Darwin" ]; then \
			M=$$(stat -f %m VERSION 2>/dev/null || echo 0); \
			if [ "$$M" -gt 0 ] 2>/dev/null; then D=$$(date -r "$$M" +%F); else D=""; fi; \
		else \
			ME=$$(stat -c %Y VERSION 2>/dev/null || echo 0); \
			if [ "$$ME" -gt 0 ] 2>/dev/null; then D=$$(date -d "@$$ME" +%F); else D=""; fi; \
		fi; \
		if [ -z "$$D" ]; then \
			D=$$(git log -1 --date=short --format=%cd -- VERSION 2>/dev/null); \
		fi; \
		if [ -z "$$D" ]; then D=$$(date +%F); fi; \
		echo "$$D" > RELEASE_DATE; \
		echo "[gen] RELEASE_DATE = $$D (from VERSION last modified)"; \
	}
	@cp RELEASE_DATE $(IBP_DIR)/IBP_RELEASE_DATE
	@cp RELEASE_DATE $(ITP_DIR)/ITBP_RELEASE_DATE
	@echo "Initialized per-algorithm metadata under $(ALG_DIR)/{iBP,iTBP}/"

help:
	@echo "Targets: build (default) | all | libs | clean | version | init-algo-meta | print-<VAR>"

# ---------------- Auto dependencies ----------------
-include $(DEPS)

