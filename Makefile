
all:  main.pdf

%.pdf : %.rmd $(wildcard *.rmd) _output.yml
	Rscript -e "system.time(rmarkdown::render('$<'))"

%.docx : %.rmd $(wildcard *.rmd) 
	
view: 
	open -a Skim main.pdf 
