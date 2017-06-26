CXX=$(shell root-config --cxx --cflags) -Wall
LIBS=$(shell root-config --libs)

%.o : %.cc
	$(CXX) -c $<

all: runPrepareEvents runMassFit runBoostAngles runReshuffleNch costh_mass

runPrepareEvents: runPrepareEvents.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runMassFit: runMassFit.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runBoostAngles: runBoostAngles.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runReshuffleNch: runReshuffleNch.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

costh_mass: massFits_costh.C
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore


clean:
	rm runPrepareEvents runMassFit runBoostAngles runReshuffleNch*.o
