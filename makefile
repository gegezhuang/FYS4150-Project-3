# https://gist.githubusercontent.com/mauriciopoppe/de8908f67923091982c8c8136a063ea6/raw/2a8a1c2c849bf8d931cc1aa4c088290c4037a02f/~Makefile

# ifndef VERBOSE
# .SILENT:
# endif

### List of all plots and data used directly in project.tex ###
ALL_TARGET_FILES = plots/some_plot.pdf

### Files used for compiling C++ or latex ###
HEADERS := $(wildcard include/*.hpp)
SOURCES := $(wildcard src/*.cpp)
TEX_FILES = $(wildcard latex/*.tex)

### C++ compiler and python runner ###
CXX ?= g++
PYTHON ?= python3.9

### C++ flags ###
GENERAL_FLAGS = -std=c++11 -larmadillo
DEBUG_FLAGS = -Wall -Wextra -g
INCLUDES = -I include

### Compiling C++ files ###
debug: $(HEADERS) $(SOURCES)
	$(CXX) $(SOURCES) $(INCLUDES) -o debug $(GENERAL_FLAGS) $(DEBUG_FLAGS)
	gdb --args ./debug

main: $(HEADERS) $(SOURCES)
	$(CXX) $(SOURCES) $(INCLUDES) -o main $(GENERAL_FLAGS)

test: test.cpp $(HEADERS) $(SOURCES)
	$(CXX) test.cpp $(SOURCES) $(INCLUDES) -o test $(GENERAL_FLAGS)
	./test

### Creating folders ###
plots:
	mkdir plots

data:
	mkdir data

### Managing virtual python environment ###
venv:
	$(PYTHON) -m venv venv
	(source venv/bin/activate)
	pip install -r requirements.txt

.PHONY: activate_venv
activate_venv: venv
	(source venv/bin/activate)

### Compile pdf from LaTeX ###
project.pdf: $(ALL_TARGET_FILES) $(TEX_FILES)
	cd latex && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
	mv latex/main.pdf project.pdf

preview: $(TEX_FILES)
	cd latex && latexmk -pdf -pvc main.tex

### General commands ###
.PHONY: all
all: clean project.pdf

.PHONY: clean
clean:
	@echo "Deleting compiled files, as well as the data, plots and venv folders"
	@rm -f main
	@rm -f test
	@rm -f debug
	@rm -f latex/main.pdf
	@rm -f -rf data
	@rm -f -rf plots
	@rm -f -rf venv
	@rm -f latex/main.log
	@rm -f latex/main.aux
	@rm -f latex/texput.log
	@rm -f latex/main.out
	@rm -f latex/main.bbl
	@rm -f latex/main.blg
	@rm -f latex/mainNotes.bib
	@rm -f latex/amsmath.aux
	@rm -f latex/main.fdb_latexmk
	@rm -f latex/main.fls

### Data ###
data/some_data.csv: data main
	touch data/some_data.csv

### Plots ###
plots/some_plot.pdf: activate_venv plots data/some_data.csv
	touch plots/some_plot.pdf
