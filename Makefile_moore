CC = g++

CFLAGS = -std=c++0x -g -Wall -Wextra

LIBS = -I /home/tools/fsl/extras/include/boost/
LIBS += -lboost_regex
#LIBS += -L /usr/include/boost/stage/lib
#LIBS += -l libboost_system


sa : main.o readFiles.o readScl.o
	$(CC) $(CFLAGS)  $(LIBS) -o sa main.o readFiles.o  readScl.o

main.o : main.h readFiles.h readScl.h main.cpp
	$(CC) $(CFLAGS) $(LIBS) -c main.cpp

readFiles.o : readFiles.h readScl.h readFiles.cpp
	$(CC) $(CFLAGS) $(LIBS) -c readFiles.cpp 

#readNets.o : readFiles.h readNets.h readNets.cpp
#	$(CC) $(CFLAGS) -c readNets.cpp

readScl.o : readFiles.h readScl.h readScl.cpp
	$(CC) $(CFLAGS) $(LIBS) -c readScl.cpp 

clean :
	-rm *.o
