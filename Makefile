# Root Makefile for pysic-rs

.PHONY: help test build build-release docs clean

help:
	@echo "pysic-rs Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  test          Run all Rust tests"
	@echo "  build         Development build (maturin develop)"
	@echo "  build-release Release build (maturin build --release)"
	@echo "  docs          Build Sphinx documentation"
	@echo "  docs-live     Live-reload documentation server"
	@echo "  clean         Remove build artifacts"
	@echo ""

test:
	cargo test

build:
	maturin develop

build-release:
	maturin develop --release

docs:
	cd docs && $(MAKE) html

docs-live:
	cd docs && $(MAKE) livehtml

clean:
	cargo clean
	rm -rf docs/build
	rm -rf target/wheels
	rm -rf *.egg-info
