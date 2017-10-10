include ~/phys_utils/make_setup_heplx.mk

all: runPolPPD

runPolPPD:
	g++ $(shell $(ROOT_CONFIG_BIN) --libs --cflags) -O2 -I$(HOME)/phys_utils/ -o runPolPPD runPolPPD.cc
