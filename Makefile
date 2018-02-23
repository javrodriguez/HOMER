SUBDIR = cpp

all:
	$(MAKE) -C $(SUBDIR)

clean:
	$(MAKE) -C $(SUBDIR) clean


.PHONY: all $(SUBDIRS)
