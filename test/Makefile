ODIR       = obj

TTIDIR     = $(CMSSW_BASE)/src/TopTagger/TopTagger/include
TPIDIR     = $(CMSSW_BASE)/src/TopTagger/CfgParser/include

TTSDIR     = $(CMSSW_BASE)/src/TopTagger/TopTagger/src
TPSDIR     = $(CMSSW_BASE)/src/TopTagger/CfgParser/src

SDIR = $(CMSSW_BASE)/src/SusyAnaTools/Tools
IDIR = $(CMSSW_BASE)/src/TopTagTest/src

CXX        = g++

CXXFLAGS  += -I$(CMSSW_BASE)/src -std=c++0x -I$(CMSSW_RELEASE_BASE)/src -I$(shell scram tool info fastjet | grep '^INCLUDE' | sed 's/INCLUDE=//')
## Optimization flag
CXXFLAGS += -g
## Enable the maximun warning
## Optimization/debug flag
CXXFLAGS += -g #-O3
## Enable for maximum warning
#CXXFLAGS += -Wall -Wextra -Weffc++

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

LIBS       = $(shell root-config --glibs)

MT2LIB = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2
TOPTAGLIB = -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lTopTaggerTopTagger
CFGPARSELIB =  -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lTopTaggerCfgParser 

PROGRAMS = remnantSys

LIBRARIES = TopTagger

all: mkobj sampPyWrap $(PROGRAMS)


mkobj:
	@mkdir -p obj

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(SDIR)/samples.h 
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(TTSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TTSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TPSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(TPSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -fPIC -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

remnantSys: $(ODIR)/remnantSys.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o  $(ODIR)/SATException.o
	$(LD) $^ $(LIBS) $(MT2LIB) $(TOPTAGLIB) $(CFGPARSELIB) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) $(LIBRARIES) core 

-include $(ODIR)/*.d