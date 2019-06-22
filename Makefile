CC := g++
CFLAGS := -g -std=c++17 -O3 -pthread -Wall -Wextra
LIBS := -lboost_regex

EAGLEFILE := test1.brd
BOOKSHELFPRFX := testing1

sa : main.o readFiles.o
	$(CC) $(CFLAGS) $(LIBS) -o sa main.o readFiles.o  #readScl.o

main.o : main.h readFiles.h main.cpp
	$(CC) $(CFLAGS) -c main.cpp

readFiles.o : readFiles.h readFiles.cpp
	$(CC) $(CFLAGS) -c readFiles.cpp

#readScl.o : readFiles.h readScl.h readScl.cpp
#	$(CC) $(CFLAGS) -c readScl.cpp

graphs:
	echo creating graphs

plot:
	echo creating plots
	#circuitname = 'bm3'
	#plfile = 'bm3.pl'
	#outfile = 'tst'
	echo $(circuitname)
	python3 make_plots.py --brd bm3 --pl bm3.pl --out tst --reports reports

animation_plot:
	circuitname = bm3
	figname = plot.gif
	echo creating animation
	python3 make_plots.py

animation_graphs:
	circuitname=bm3
	figname=plot.gif
	echo creating animation
	python3 make_plots.py

convert:
	if [ ! -d "./eagle2bookshelf" ]; then \
		echo downloading eagle2bookshelf tool; \
		git clone https://github.com/djmerrill/eagle2bookshelf.git; \
	fi
	if [ ! -f "./eagle2bookshelf/DRU.py" ]; then \
			echo downloading DRU library; \
			wget https://raw.githubusercontent.com/NVSL/Swoop/master/Swoop/DRU.py -P ./eagle2bookshelf/; \
	fi
	if [ ! -d "./cache" ]; then \
			echo "creating cache directory"; \
			mkdir "cache"; \
	fi
	if [ ! -d "./cache/designs" ]; then \
			echo "creating cache/designs directory"; \
			mkdir "./cache/designs"; \
	fi
	echo "converting from eagle to bookshelf"

	python eagle2bookshelf/eagle2bookshelf2012.py --brd $(EAGLEFILE) --output_prfx $(BOOKSHELFPRFX) --userid annealer
	mv $(BOOKSHELFPRFX).nodes ./cache/designs
	mv $(BOOKSHELFPRFX).nets ./cache/designs
	mv $(BOOKSHELFPRFX).pl ./cache/designs
	mv $(BOOKSHELFPRFX).wts ./cache/designs


install:
	if [ ! -d "./cache" ]; then \
			echo "creating cache directory"; \
			mkdir "cache"; \
	fi
	if [ ! -d "./cache/img" ]; then \
			echo "creating cache/img directory"; \
			mkdir "./cache/img"; \
	fi
	if [ ! -d "./cache/designs" ]; then \
			echo "creating cache/designs directory"; \
			mkdir "./cache/designs"; \
	fi
	if [ ! -d "./cache_rudy" ]; then \
			echo "creating rudy cache directory"; \
			mkdir "./cache_rudy"; \
	fi
	if [ ! -d "./reports" ]; then \
			echo "creating reports directory"; \
			mkdir "reports"; \
	fi

clean :
	-rm *.o
	-rm cache/*.pl
	-rm cache/designs/*.pl
	-rm cache/img/*.png
	-rm cache_rudy/*.txt
	-rm reports/*.txt
	-rm wl.txt
	-rm oa.txt
	-rm cost.txt
	-rm anim.gif
