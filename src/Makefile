CC := g++
CFLAGS := -g -std=c++17 -O3 -pthread # -Wall -Wextra
LIBS := -lboost_regex

EAGLEFILE := test1.brd
BOOKSHELFPRFX := testing1

annealer : main.o readFiles.o
	$(CC) $(CFLAGS) $(LIBS) -o annealer main.o readFiles.o  #readScl.o

main.o : main.hpp readFiles.hpp main.cpp
	$(CC) $(CFLAGS) -c main.cpp

readFiles.o : readFiles.hpp readFiles.cpp
	$(CC) $(CFLAGS) -c readFiles.cpp

#readScl.o : readFiles.hpp readScl.hpp readScl.cpp
#	$(CC) $(CFLAGS) -c readScl.cpp
