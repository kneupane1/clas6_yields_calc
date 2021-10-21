ROOTLIBS	= $(shell root-config --libs)
CXX = $(shell root-config --cxx)
CXXFLAGS = -O3 -march=native -fPIC -w -g $(shell root-config --cflags)

SRCDIR   = src
OBJDIR   = obj
LIBDIR 	 = lib
BINDIR	 = bin
TARGET 	 = Event

SRC = $(wildcard $(SRCDIR)/*.cpp)

OBJECTS  = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
LIBOUT =	$(LIBDIR)/libEvent.so

.PHONY: all clean

all: obj_dir $(BINDIR)/$(TARGET)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

obj_dir:
	@mkdir -p $(OBJDIR)

lib_dir:
	@mkdir -p $(LIBDIR)

bin_dir:
	@mkdir -p $(BINDIR)

lib:	obj_dir lib_dir $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared $(OBJECTS) $(ROOTLIBS) -o $(LIBOUT)

$(MAIN): $(OBJDIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BINDIR)/$(TARGET): bin_dir $(OBJECTS) $(MAIN)
	$(CXX) $(OBJECTS) $(MAIN) $(CXXFLAGS) $(ROOTLIBS) -o $(BINDIR)/$(TARGET)

clean:
	-rm -f $(OBJECTS) $(LIBOUT)
