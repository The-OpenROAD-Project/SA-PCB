all:
	+$(MAKE) -C src
	mv src/annealer ./bin

graphs:
	echo creating graphs

circuitname := './designs/bm1'
plfile := './bin/cache/10.pl'
outfile := './tst'

plot:
	echo creating plots
	python3 src/py_utils/make_plots.py --brd $(circuitname)  --pl $(plfile) --out $(outfile) --reports reports

animate_plot:
	circuitname = bm3
	figname = plot.gif
	echo creating animation
	python3 src/py_utils/make_plots.py

animate_graphs:
	circuitname=bm3
	figname=plot.gif
	echo creating animation
	python3 src/py_utils/make_plots.py

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
	cp -r cpp-taskflow/taskflow .
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
	-rm src/*.o
	-rm cache/*.pl
	-rm cache/designs/*.pl
	-rm cache/img/*.png
	-rm cache_rudy/*.txt
	-rm reports/*.txt
