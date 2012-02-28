CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

SOURCE =	src

OUTPUT =	bin

LIBS =

all:	SINAH decompose

SINAH:	$(SOURCE)/SINAH.o
	$(CXX) -o $(OUTPUT)/SINAH $(SOURCE)/SINAH.o $(LIBS)
	
decompose:	$(SOURCE)/decompose.o
	$(CXX) -o $(OUTPUT)/decompose $(SOURCE)/decompose.o $(LIBS)



clean:
	rm -f $(OBJS) $(TARGET)
