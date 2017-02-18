HFILES=EdgeBoxGenerator.h \
extern.h \
Image.h	\
model.h \
Option.h \
ttq.h

OBJS=EdgeBoxGenerator.o \
extern.o \
Image.o \
main.o \
model.o \
Option.o \
wrappers.o

CFLAGS=-O2 -I. -I/usr/include/arm-linux-gnueabihf -I/usr/local/include -I/usr/include -mfloat-abi=hard -mfpu=neon -g3 -DDEBUG
LDFLAGS=-lpthread $(shell pkg-config --libs opencv)  #-lopencv


all: $(OBJS) $(HFILES)
	g++ $(CFLAGS) -o myedgeboxes $(OBJS) $(LDFLAGS) 

%o:%cpp $(HFILES)
	g++ $(CFLAGS) -c $<

clean:
	rm -rf *.o
	rm -rf myedgeboxes
