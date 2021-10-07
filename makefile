# https://gist.githubusercontent.com/mauriciopoppe/de8908f67923091982c8c8136a063ea6/raw/2a8a1c2c849bf8d931cc1aa4c088290c4037a02f/~Makefile

# ifndef VERBOSE
# .SILENT:
# endif

### List of all plots and data used directly in project.tex ###
ALL_TARGET_FILES = plots/some_plot.pdf

### C++ compiling variables ###
HEADERS := $(wildcard include/*.hpp)
SOURCES := $(wildcard src/*.cpp)

CXX ?= g++
PYTHON ?= python3.9

GENERAL_FLAGS = -std=c++11 -larmadillo
DEBUG_FLAGS = -Wall -Wextra -g
INCLUDES = -I include

### Compiling C++ files ###
debug: src/main.cpp $(HEADERS) $(SOURCES)
	$(CXX) src/main.cpp $(SOURCES) $(INCLUDES) -o debug $(GENERAL_FLAGS) $(DEBUG_FLAGS)
	gdb --args ./debug

main: src/main.cpp $(HEADERS) $(SOURCES)
	$(CXX) src/main.cpp $(SOURCES) $(INCLUDES) -o main $(GENERAL_FLAGS)

test: src/test.cpp $(HEADERS) $(SOURCES)
	$(CXX) src/test.cpp $(SOURCES) $(INCLUDES) -o test $(GENERAL_FLAGS)
	./test

### Creating folders ###
plots:
	mkdir plots

data:
	mkdir data

### Managing virtual python environment ###
venv:
	$(python) -m venv venv
	(source venv/bin/activate)
	pip install -r requirements.txt

.PHONY activate_venv
activate_venv: venv
	(source venv/bin/activate)

### Compile pdf from LaTeX ###
project.pdf: $(all_target_files)
	cd latex && pdflatex project.tex && pdflatex project.tex
	mv latex/project.pdf project.pdf

### General commands ###
.PHONY: all
all: clean project.pdf

.PHONY: clean
clean:
	@echo "Deleting compiled files, as well as the data, plots and venv folders"
	@rm -f main
	@rm -f test
	@rm -f debug
	@rm -f -rf data
	@rm -f -rf plots
	@rm -f -rf venv
	@rm -f latex/project.log
	@rm -f latex/project.aux
	@rm -f latex/texput.log
	@rm -f latex/project.out
	@rm -f latex/projectNotes.bib
	@rm -f latex/amsmath.aux

### Data ###
data/some_data.csv: data main
	./main some_data

### Plots ###
plots/some_plot.pdf: activate_venv plots data/some_data.csv
	python src/plot_some_plot.py
