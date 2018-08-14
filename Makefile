clean:
	rm -v src/*.so src/*.o
	rm -v R/RcppExports.R
	rm -v src/RcppExports.cpp

build:
	Rscript .compileAttributes.R

install:
	R CMD INSTALL ../spectralGraphTopology
