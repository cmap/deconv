opts = --standalone --template=templates/boyd.latex \
      --natbib --bibliography=dpeak.bib \
      --from markdown -t latex+raw_tex

md_files = main.md 01-intro.md 02-meth.md 03-results.md 04-discuss.md 05-acknow.md 10-figs.md 20-support.md

all: draft.pdf view

draft.md : $(md_files)
	cat $^ > draft.md

draft.pdf: draft.md
	pandoc $(opts) -o draft.tex draft.md
	pdflatex draft.tex
	bibtex draft
	pdflatex draft.tex
	pdflatex draft.tex
	wc -w 0*.md 10-figs.md

view:
	open -a skim draft.pdf

clean:
	rm draft.*
