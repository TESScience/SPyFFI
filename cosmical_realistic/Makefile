.PHONY: all clean

all: _cosmical.so

_cosmical.so:
	python setup.py build_ext --inplace

clean:
	rm -rf build/ *.so
