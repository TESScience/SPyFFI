TESTDIR := $(shell pwd)
# where the Python scripts can be found (i.e. the directory containing SPyFFI/ and zachopy/)
export PYTHONPATH:=$(TESTDIR)/../..

# the root directory of SPyFFI (needed for some file-finding by the code)
export SPYFFIPATH=$(TESTDIR)/..

# the directory where you want to store TESS data (inputs/intermediates/outputs)
export SPYFFIDATA=$(TESTDIR)/SPyFFIdata/

.PHONY: tests smallone clean

tests: smallone

smallone: $(SPYFFIDATA)/inputs
	./smallone.py

$(SPYFFIDATA):
	mkdir $(SPYFFIDATA)

$(SPYFFIDATA)/inputs: $(SPYFFIDATA)
	(cd $(SPYFFIDATA); \
		wget https://www.dropbox.com/s/0e4c2uk34phv4qx/SPyFFI_coreinputs.tar.gz; \
		tar xzf SPyFFI_coreinputs.tar.gz )

clean:
	rm -rf $(SPYFFIDATA)
