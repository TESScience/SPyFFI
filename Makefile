all:

clean:
	rm -f $(shell find . -name "*.pyc")
	rm -rf build/ dist/ SPyFFI.egg-info/
