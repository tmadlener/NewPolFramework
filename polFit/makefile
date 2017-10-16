include ~/phys_utils/make_setup_heplx.mk

all: runPolPPD runPolGen runPolSub

runPolPPD:
	g++ $(shell $(ROOT_CONFIG_BIN) --libs --cflags) -O2 -I$(HOME)/phys_utils/ -o runPolPPD runPolPPD.cc

runPolGen:
	g++ $(shell $(ROOT_CONFIG_BIN) --libs --cflags) -O2 -I$(HOME)/phys_utils/ -o runPolGen runPolGen.cc

runPolSub:
	g++ $(shell $(ROOT_CONFIG_BIN) --libs --cflags) -O2 -I$(HOME)/phys_utils/ -o runPolSub runPolSub.cc
