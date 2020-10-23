r_scripts = $(wildcard 0*.R)
all: $(r_scripts:.R=.html)

%.Rmd : %.R
	Rscript -e "knitr::spin('$<', knit = F)" 

%.html : %.Rmd
	Rscript -e "rmarkdown::render('$<')"

clean:
	rm *.html tab*.md

view:
	open *.html 
