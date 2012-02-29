CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

SOURCE =	src

OUTPUT =	bin

LIB_DIRS =

LIBS = 		hdf5 hdf5_hl OGDF pthread

LDFLAGS += $(foreach librarydir,$(LIB_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))


all:	SINAH decompose

SINAH:	$(SOURCE)/SINAH.o
	$(CXX) -o $(OUTPUT)/SINAH $(SOURCE)/SINAH.o $(LDFLAGS)
	
decompose:	$(SOURCE)/decompose.o $(SOURCE)/BundleGraph.o $(SOURCE)/types.o
	$(CXX) -o $(OUTPUT)/decompose $(SOURCE)/decompose.o $(SOURCE)/BundleGraph.o $(SOURCE)/types.o $(LDFLAGS)


clean:
	rm -f $(OBJS) $(TARGET)
