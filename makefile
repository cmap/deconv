
all:  main.docx \
      main.pdf \
      submission/nat_meth/cover.pdf

%.pdf : %.rmd $(wildcard *.rmd) _output.yml
	Rscript -e "system.time(rmarkdown::render('$<'))"

%.docx : %.rmd $(wildcard *.rmd) 
	Rscript -e "system.time(rmarkdown::render('$<', output_format='word_document'))"

submission/nat_meth/%.pdf: submission/nat_meth/%.rmd
	Rscript -e "system.time(rmarkdown::render('$<'))"

figures.rout : scripts/figures.R
	nohup Rscript "$<" &> "$@" 

	
view: 
	open -a Skim main.pdf 
