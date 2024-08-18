CXX = g++
CXXFLAGS = -O3 -Wall

all : test

test: bin/test

clean :
	- rm -r bin
	- rm -r obj

obj :
	mkdir obj

bin :
	mkdir bin

# Compile objects
obj/write.o: writer.cpp writer.h | obj
	$(CXX) -c writer.cpp -o obj/writer.o $(CXXFLAGS)

bin/test: main.cpp obj/write.o | bin
	$(CXX) main.cpp obj/writer.o -o bin/test $(CXXFLAGS)
