SHELL = /bin/sh
#

#PRODDIR and SUBDIR are passed through MFLAGS mechanism.

all:
	@echo all target does nothing - make complete

install :
	@if [ "$(PRODDIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi 
	cp -p *.dat Makefile $(PRODDIR)/$(SUBDIR)
	for datadir in dsc decam despec filter lsst \
		prime rover sdss snap; do \
	   (mkdir $(PRODDIR)/$(SUBDIR)/$$datadir; cd $$datadir; \
	    cp *.lns *.len *.dat *.tcl $(PRODDIR)/$(SUBDIR)/$$datadir; \
	    cd ..;) \
	done;
	chmod -R 755 $(PRODDIR)/$(SUBDIR)/*

clean :
	rm -f *.bak *~ *.tkbak '#'*'#'
