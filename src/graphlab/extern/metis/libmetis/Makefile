include ../Makefile.in

CURBUILDDIR = $(LIBBUILDDIR)

OBJS      = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(wildcard *.c))
HEADERS   = $(wildcard *.h) $(TOPDIR)/include/metis.h  $(wildcard $(GKLIBINCDIR)/*.h)

GKLIBOBJS = $(wildcard $(GKLIBBUILDDIR)/*$(OBJEXT))

$(BUILDDIR)/libmetis$(LIBEXT): $(OBJS) $(GKLIBOBJS)
	$(AR) $(OBJS) $(GKLIBOBJS)
	$(RANLIB) 

clean:
	rm -f $(OBJS)

realclean:
	rm -f $(OBJS) ; rm -f $(BUILDDIR)/libmetis$(LIBEXT)

dist:
	mkdir $(TOPDIR)/../$(PKGNAME) ;\
	cp -r $(TOPDIR)/* $(TOPDIR)/../$(PKGNAME) ;\
	tar -C $(TOPDIR)/../ --exclude=.svn -czf $(TOPDIR)/$(PKGNAME).tar.gz $(PKGNAME);\
	rm -rf $(TOPDIR)/../$(PKGNAME) 

$(OBJS) : $(HEADERS) ../Makefile.in Makefile


$(CURBUILDDIR)/%$(OBJEXT) : %.c
	$(CC) $(CFLAGS) $(SOURCEFILE) $(OUTPUTFILE)
