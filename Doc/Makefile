#VERSION = svn
HELPDOC = ../../dev-tools/helpdoc #-version $(VERSION)

#LATEX   = pdflatex
#LATEX2HTML = latex2html

#PDFS = user_guide.pdf
#AUXS = $(PDFS:.pdf=.aux)
#LOGS = $(PDFS:.pdf=.log)
#OUTS = $(PDFS:.pdf=.out)
#TOCS = $(PDFS:.pdf=.toc)

#default: all

#doc:  all
#all:  pdf html defs
#all:  defs
#pdf: $(PDFS)
#html: user_guide
all: input_xx.xsl INPUT_WSTAT.html INPUT_WSTAT.txt INPUT_WFREQ.html INPUT_WFREQ.txt INPUT_WESTPP.html INPUT_WESTPP.txt 

#$(PDFS): %.pdf: %.tex
# $(LATEX)  $<
# $(LATEX)  $<

clean:
	- rm -rf input_xx.xsl
	- rm -f INPUT_WSTAT.html INPUT_WSTAT.txt INPUT_WSTAT.xml
	- rm -f INPUT_WFREQ.html INPUT_WFREQ.txt INPUT_WFREQ.xml
	- rm -f INPUT_WESTPP.html INPUT_WESTPP.txt INPUT_WESTPP.xml
# - rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~


input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1; \
	fi) ; fi)

INPUT_WSTAT.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_WSTAT.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_WFREQ.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_WFREQ.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_WESTPP.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_WESTPP.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<
