# Paths
SUBMODULE_DIR = external/wfa2
PYTHON = python3
PIP = pip

.PHONY: all wfa2 install clean clean-install submodules dev-install hilbert-install

all: wfa2 install

clean-install: clean install

submodules: wfa2

# rules for installing submodules
# -------------------------------

wfa2:
	@if [ ! -f external/wfa2/Makefile ]; then \
		echo "\nWFA2 module missing. Initializing..."; \
		git submodule update --init --recursive external/wfa2; \
	fi
	@echo "\nBuilding C++ wfa2 submodule..."
	$(MAKE) -C external/wfa2 clean all

# installation rules
# ------------------

install:
	@echo "\nCompiling Cython and installing giggles via pip..."
	$(PIP) install .

# assumes that user is already in the giggles-dev environment
dev-install: submodules
	@echo "\nCompiling Cython and installing developmental giggles via pip..."
	pip install -e .[dev]

hilbert-install:
	@echo "\nSetting up environment, compiling submodules, and installing via pip..."
	module load gcc/15.1.0
	$(MAKE) submodules
	PIP_CONFIG_FILE=/software/python/pip.conf pip install -e .[dev]

# cleanup
# -------

clean:
	@echo "\nCleaning Python/Cython build artifacts..."
	rm -rf build/ dist/ *.egg-info/
	find . -name "*.pyc" -delete
	find . -name "*.so" -delete
	find . -name "__pycache__" -delete