CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

SOURCE =	src

OUTPUT =	bin

LIB_DIRS =	/usr/lib /usr/local/lib

LIBS = 		hdf5 hdf5_hl OGDF

LDFLAGS += $(foreach librarydir,$(LIB_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

#INC_DIRS = 	/opt/OGDF/_release/

#INCS = 		OGDF

#INCFLAGS += $(foreach librarydir,$(INC_DIRS),-I$(librarydir))
#INCFLAGS += $(foreach library,$(INCS),-include$(library))


all:	SINAH decompose

SINAH:	$(SOURCE)/SINAH.o
	$(CXX) -o $(OUTPUT)/SINAH $(SOURCE)/SINAH.o $(LDFLAGS)
	
decompose:	$(SOURCE)/decompose.o
	$(CXX) -o $(OUTPUT)/decompose $(SOURCE)/decompose.o $(LDFLAGS)


clean:
	rm -f $(OBJS) $(TARGET)
