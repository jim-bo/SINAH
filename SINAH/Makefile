CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		SINAH.o

LIBS =

TARGET =	SINAH

SOURCE =	src

OUTPUT =	bin

$(TARGET):	$(SOURCE)/$(OBJS)
	$(CXX) -o $(OUTPUT)/$(TARGET) $(SOURCE)/$(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
