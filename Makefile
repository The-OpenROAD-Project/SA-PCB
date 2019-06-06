CC := g++
CFLAGS := -g -std=c++17 -O3 -pthread -Wall -Wextra
LIBS := -lboost_regex

sa : main.o readFiles.o readScl.o
	$(CC) $(CFLAGS) $(LIBS) -o sa main.o readFiles.o  readScl.o

main.o : main.h readFiles.h readScl.h main.cpp
	$(CC) $(CFLAGS) -c main.cpp

readFiles.o : readFiles.h readScl.h readFiles.cpp
	$(CC) $(CFLAGS) -c readFiles.cpp

#readNets.o : readFiles.h readNets.h readNets.cpp
#	$(CC) $(CFLAGS) -c readNets.cpp

readScl.o : readFiles.h readScl.h readScl.cpp
	$(CC) $(CFLAGS) -c readScl.cpp

install :
	if [ ! -d "./cache" ]; then \
			echo "creating cache directory"; \
			mkdir "cache"; \
	fi
	if [ ! -d "./cache_rudy" ]; then \
			echo "creating rudy cache directory"; \
			mkdir "cache_rudy"; \
	fi
	if [ ! -d "./reports" ]; then \
			echo "creating reports directory"; \
			mkdir "reports"; \
	fi

clean :
	-rm *.o
	-rm *.pl
	-rm cache/*.pl
	-rm cache/img/*.png
	-rm cache_rudy/*.txt
	-rm reports/*.txt
