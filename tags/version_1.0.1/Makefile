R_BIN=R
FILE_NAME=peakROTS_1.0.1

all:
	$(R_BIN) CMD build peakROTS
	$(R_BIN) CMD Rd2pdf --no-preview --pdf peakROTS
	mv peakROTS.pdf ${FILE_NAME}.pdf

install:
	$(R_BIN) CMD INSTALL ${FILE_NAME}.tar.gz
	@echo
	@echo "Run $(R_BIN), then test:"
	@echo "library(peakROTS)"
	@echo "help(peakROTS)"


clean:
	rm -f ${FILE_NAME}.tar.gz
	rm -f ${FILE_NAME}.pdf
