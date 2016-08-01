.PHONY: all clean

all: cosmical.so

cosmical.so:
	python compile_module build_ext --inplace

clean:
	rm -rf build/ *.so
