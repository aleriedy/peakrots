R_BIN=/v/linux26_x86_64/appl/stat/R/R-2.10.1/bin/R
OLD_R_BIN=R
R_LIB_DIR=/wrk/akallio/R_libs
FILE_NAME=peakROTS_1.0.0

all:
	$(OLD_R_BIN) CMD build peakROTS
	$(R_BIN) CMD Rd2pdf --no-preview --pdf peakROTS
	mv peakROTS.pdf ${FILE_NAME}.pdf

install:

	mkdir -p $(R_LIB_DIR)
	$(OLD_R_BIN) CMD INSTALL -l $(R_LIB_DIR) ${FILE_NAME}.tar.gz
	@echo
	@echo "Run export R_LIBS_USER=$(R_LIB_DIR) and $(OLD_R_BIN), then test:"
	@echo "library(peakROTS)"
	@echo "help(peakROTS)"


clean:
	rm -f ${FILE_NAME}.tar.gz
	rm -f ${FILE_NAME}.pdf
