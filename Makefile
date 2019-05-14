
all:  main.pdf

%.pdf : %.rmd $(wildcard *.rmd) _output.yml
	Rscript -e "system.time(rmarkdown::render('$<'))"
	
view: 
	open -a Skim main.pdf 
