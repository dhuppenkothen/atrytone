# Put the directory above 'DNest4' into an environment variable
# called DNEST4_PATH for this to work.
#XSDIR =  /Users/danielahuppenkothen/work/sw/heasoft-6.17/Xspec/x86_64-apple-darwin14.5.0
XSDIR = /Users/danielahuppenkothen/work/sw/heasoft-6.21/x86_64-apple-darwin16.6.0

CXX=g++
#CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
CXXFLAGS = -std=c++11 -g -O0 -Wall -Wextra -pedantic -rpath $(XSDIR)/lib
#LIBS = -ldnest4 -lpthread -lgsl -lgslcblas -L/usr/local/lib/ -lCCfits -L/opt/local/lib/ -L$(XSDIR)/lib -lcfitsio -lXSFunctions -lXSModel -lXSUtil -lXS
LIBS = -ldnest4 -lpthread -lgsl -lgslcblas -L$(XSDIR)/lib -lcfitsio -lCCfits_2.5 -lXSFunctions -lXSModel -lXSUtil -lXS  -L/opt/local/lib/ 
 
default:
#	g++ -I$(DNEST4_PATH) -I$(XSDIR)/include -I/usr/local/include/ccfits -r/usr/include/CCfits -I/opt/local/include $(CXXFLAGS) -c *.cpp
	$(CXX) $(CXXFLAGS) -I$(DNEST4_PATH) -I$(XSDIR)/include -I$(XSDIR)/include/CCfits -I/opt/local/include -c *.cpp

	$(CXX) $(CXXFLAGS) -L$(DNEST4_PATH)/DNest4/code -o main *.o $(LIBS)
	rm *.o
