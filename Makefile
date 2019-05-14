
all:  main.pdf

%.pdf : %.rmd $(wildcard *.rmd) _output.yml
	Rscript -e "system.time(rmarkdown::render('$<'))"

%.docx : %.rmd $(wildcard *.rmd) 
	Rscript -e "system.time(rmarkdown::render('$<', output_format='word_document'))"
	
view: 
	open -a Skim main.pdf 
