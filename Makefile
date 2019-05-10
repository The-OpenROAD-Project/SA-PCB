CC = g++
CFLAGS = -g -Wall -Wextra

sa : main2.o readFiles.o readScl.o
	$(CC) $(CFLAGS) -o sa main2.o readFiles.o  readScl.o

main2.o : main.h readFiles.h readScl.h main2.cpp
	$(CC) $(CFLAGS) -c main2.cpp

readFiles.o : readFiles.h readScl.h readFiles.cpp
	$(CC) $(CFLAGS) -c readFiles.cpp

#readNets.o : readFiles.h readNets.h readNets.cpp
#	$(CC) $(CFLAGS) -c readNets.cpp

readScl.o : readFiles.h readScl.h readScl.cpp
	$(CC) $(CFLAGS) -c readScl.cpp

clean :
	-rm *.o
