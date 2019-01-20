TMPDIR := $(shell mktemp -d)
PKGDIR := "${HOME}/Development/stochvol"

## R HOME
R_HOME := $(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 		$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := 		$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
#RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 		$(shell echo 'library(Rcpp, lib.loc="~/R/dev_version/"); Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'library(Rcpp, lib.loc="~/R/dev_version/"); Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'library(RInside, lib.loc="~/R/dev_version/"); RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'library(RInside, lib.loc="~/R/dev_version/"); RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## RcppArmadillo headers
RCPPARMAINCL :=		$(shell echo 'library(RcppArmadillo, lib.loc="~/R/dev_version/"); RcppArmadillo:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)	

## compiler etc settings used in default make rules
CXX := 			$(shell $(R_HOME)/bin/R CMD config CXX)
CPPFLAGS := 		-Wall $(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS := 		$(RCPPFLAGS) $(RCPPINCL) $(RCPPARMAINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)
LDLIBS := $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

sources := 		$(wildcard src/*.cpp)
objects := $(sources:.cpp=.o)

install:
	#R --vanilla -e "devtools::document(\"${PKGDIR}\")"
	R CMD INSTALL --no-multiarch --with-keep.source --library=${HOME}/R/under_development ${PKGDIR}

debug: src/stochvol.so

%.o: %.cpp
	$(CXX) -xc++ -c -fPIC -std=c++11 $(CPPFLAGS) $(CXXFLAGS) -g -fdebug-info-for-profiling -O0 -o $@ $<

src/stochvol.so: $(objects)
	$(CXX) -shared -std=c++11 -g -fdebug-info-for-profiling $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS) -L./src/ -Wl,-rpath,/home/dhosszejni/Development/stochvol/src/ -O0 $^ -o $@

build:
	R --vanilla -e "devtools::document(\"${PKGDIR}\")"
	'/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore --quiet  \
		CMD INSTALL ${PKGDIR}  \
		--library="${TMPDIR}" --no-R --no-data  \
		--no-help --no-demo --no-inst --no-docs --no-exec --no-multiarch  \
		--no-test-load --preclean

clean:
	rm -fv src/*.o src/*.so

test:
	$(CXX) -xc++ -std=c++11 $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS) -g -fprofile-instr-generate -fdebug-info-for-profiling -O0 -o test_cpp test.cpp $(sources)
