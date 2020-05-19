MANUSCRIPT	= alife2020_sustainability
TEXFILES	= $(MANUSCRIPT:%=%.tex)
PSFILES		= $(MANUSCRIPT:%=%.ps)
PDFFILES	= $(MANUSCRIPT:%=%.pdf)
PDFLABSHEETS	= $(LABSHEETS:%=%.pdf)

all : $(PDFFILES)

$(TEXFILES) : alifeconf.sty al20sustainfigs.Rout

clean :
	rm -f *~ *.dvi *.aux *.out *.log *.bbl *.blg $(PSFILES) $(PDFFILES)

.PHONY : testslide allslides allsheets all clean tgz
.PRECIOUS : %.dvi %.ps %.Rout
.SECONDARY :

%.pdf : %.ps
	ps2pdf -dAutoRotatePages=/None $< $@
#	ps2pdf -dEmbedAllFonts=true -dPDFSETTINGS=/prepress -dAutoRotatePages=/None $< $@

%.ps : %.dvi
	dvips -Z -P pdf -o $@ $<
#	dvips -Z -P amz -P cmz -o $@ $<

%.eps : %.xwd
	xwdtopnm $< | pnmtops -noturn > $@

%.eps : %.jpg
	jpegtopnm $< | pnmtops -noturn > $@

%.eps : %.gif
	giftopnm $< | pnmtops -noturn > $@

%.eps : %.png
	pngtopnm -mix -background '#FFF' $< | pnmtops -noturn > $@

%.eps : %.svg
	inkscape $<  --export-eps=$@


%.Rout : %.R
	R --vanilla CMD BATCH $<


%.dvi : %.tex
	latex $< < /dev/null
	if grep -q bibdata $*.aux ; then bibtex $* ; fi
	latex $<
	latex $<
	latex $<

.SUFFIXES :
