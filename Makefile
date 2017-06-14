CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

%.o : %.cc
	$(CXX) -c $<

all: runPrepareEvents runMassFit runBoostAngles runReshuffleNch

runPrepareEvents: runPrepareEvents.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
	
runMassFit: runMassFit.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
	
runBoostAngles: runBoostAngles.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runReshuffleNch: runReshuffleNch.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
	
									
clean: 
	rm runPrepareEvents runMassFit runBoostAngles runReshuffleNch*.o