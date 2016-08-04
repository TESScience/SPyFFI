.PHONY: test clean

TEST_OUT=SPyFFIdata/outputs/18h00m00s+66d33m39s_smallone/1800s/sub400x400/simulated_18h00m00s+66d33m39s_sub400x400_000000.fits
test: clear-tests $(TEST_OUT) light_curve_tests


###################### Virtual Environment ######################

.PHONY: install pre-install uninstall reinstall fix-matplotlib-config
VIRTUAL_ENV=$(CURDIR)/venv
PYTHON=$(VIRTUAL_ENV)/bin/python2.7
PIP=$(VIRTUAL_ENV)/bin/pip2.7
SITE_PACKAGES=$(VIRTUAL_ENV)/lib/python2.7/site-packages
VERSION=$(shell python -c "$(shell grep '^VERSION' $(CURDIR)/../setup.py) ; print VERSION")
INSTALL=$(SITE_PACKAGES)/SPyFFI-$(VERSION)-py2.7.egg

$(VIRTUAL_ENV):
	virtualenv --python=python2.7 $@
	@[ -d $@ ]
	touch $@

$(SITE_PACKAGES)/numpy: $(VIRTUAL_ENV)
	$(PIP) install --quiet $(notdir $@)
	@[ -d $@ ]
	touch $@

fix-matplotlib-config: $(SITE_PACKAGES)/matplotlib
	$(eval \
	MATPLOTLIBRC:=$(shell $(PYTHON) -c \
		"import matplotlib as mpl ; \
		 import os.path ; \
		 print os.path.join(mpl.get_configdir(),'matplotlibrc')"))
	mkdir -p "$(dir $(MATPLOTLIBRC))"
	@# Hack to deal with http://stackoverflow.com/a/21789908/586893
	@if [ "$(shell uname)" = "Darwin" ] ; then \
		echo "Writing 'backend: TkAgg' to $(MATPLOTLIBRC)"; \
		echo 'backend: TkAgg' > $(MATPLOTLIBRC); \
	fi
	@# Hack to deal with http://stackoverflow.com/a/3054314/586893
	@if [ "$(shell uname)" = "Linux" ] ; then \
		echo "Writing 'backend: agg' to $(MATPLOTLIBRC)"; \
		echo 'backend: agg' > $(MATPLOTLIBRC); \
	fi

$(INSTALL): $(SITE_PACKAGES)/numpy $(VIRTUAL_ENV)
	make fix-matplotlib-config
	make pre-install
	(cd ..; $(PYTHON) setup.py install)
	[ -e $(SITE_PACKAGES)/SPyFFI-$(VERSION)-py2.7-*.egg ]
	rm -f $@
	ln -s $(SITE_PACKAGES)/SPyFFI-$(VERSION)-py2.7-*.egg $@
	[ -L $@ ]

$(SITE_PACKAGES)/%: $(SITE_PACKAGES)/numpy $(VIRTUAL_ENV)
	$(PIP) install --quiet $(notdir $@)
	@[ -d $@ ]
	touch $@

pre-install: $(SITE_PACKAGES)/numpy
	$(PIP) install -U --quiet pip setuptools
	$(PIP) install --quiet -r requirements.txt

install: $(INSTALL)

uninstall: $(INSTALL)
	echo y | $(PIP) uninstall SPyFFI


###################### Testing ######################

light_curve_tests: light_curve_draw_rotation_test light_curve_draw_transit_test

SPyFFIdata/inputs:
	mkdir -p $(dir $@)
	wget -qO- https://www.dropbox.com/s/0e4c2uk34phv4qx/SPyFFI_coreinputs.tar.gz | tar xzv -C $(dir $@)

$(TEST_OUT): SPyFFIdata/inputs $(INSTALL)
	# PYTHONPATH: Where the Python scripts can be found (i.e. the directory containing SPyFFI/ and zachopy/)
	# SPYFFIPATH: The root directory of SPyFFI (needed for some file-finding by the code)
	# SPYFFIDATA: The directory where you want to store TESS data (inputs/intermediates/outputs)
	LOG=INFO SPYFFIPATH=$(CURDIR)/.. SPYFFIDATA=$(CURDIR)/SPyFFIdata/ $(PYTHON) ./scripts/smallone.py

light_curve_%_test: light_curve_%.json
	./scripts/json_diff $< reference/$<

light_curve_%.json: ./scripts/light_curve_%.py $(INSTALL)
	LOG=INFO $(PYTHON) $< > $@

clear-tests:
	rm -rf SPyFFIdata/outputs *.out

clean:
	rm -rf SPyFFIdata $(VIRTUAL_ENV)
