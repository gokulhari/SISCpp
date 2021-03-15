DOCS=index ex01 ex02 ex03 ex04 ex05 ex06 ex07 ex08 ex09 ex10 ex11 ex12 ex13 download feedback
DOCS1=$(addprefix jemdocs/, $(DOCS))
HDOCS=$(addsuffix .html, $(DOCS))
#PHDOCS=$(addprefix html/, $(HDOCS))
PHDOCS=$(HDOCS)

.PHONY : docs
docs : $(PHDOCS)

.PHONY : update
update : $(PHDOCS)
	@echo -n 'Copying to server...'
	# insert code for copying to server here.
	@echo ' done.'

%.html : jemdocs/%.jemdoc MENU
	python2.7 jemdoc.py -c bibover.conf -o $@ $<

.PHONY : clean
clean :
	rm -f *.html
