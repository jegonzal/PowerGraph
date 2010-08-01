#*************************************************************************
# Global flags
#*************************************************************************
gdb         = yes
debug       = no
memdbg      = no
openmp      = no
pcre        = no
x86compiler = gcc



#*************************************************************************
# System-specific compilation flags
#*************************************************************************
# Get some basic information about the system that you are working on
cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname)
ifeq ($(findstring CYGWIN, $(systype)),CYGWIN)
  systype = MSWIN
#  systype = CYGWIN
  cputype = x86
endif



ifeq ($(systype),MSWIN)
  #-------------------------------------------------------------------
  # These defs are very much Visual Studio Specific
  #-------------------------------------------------------------------
  #Compiler information
  CC = cl 
  OPTFLAGS = /Ox 
  COPTIONS = -DWIN32 -DMSC -D_CRT_SECURE_NO_DEPRECATE 

  #Linker information
  LDOPTIONS = 
  LD = /cygdrive/c/Program\ Files/Microsoft\ Visual\ Studio\ 8/VC/BIN/link

  #Library creation information
  AR = lib /NODEFAULTLIB:LIBC /OUT:$@ 
  RANLIB = 

  ifeq ($(openmp),yes)
    COPTIONS  += -D__OPENMP__ /openmp 
    LDOPTIONS += /openmp
  endif

  #Compile input/output file specification
  SOURCEFILE = /c $<
  OUTPUTFILE = /Fo$@

  # Standard file extensions
  OBJEXT = .obj
  LIBEXT = .lib
  EXEEXT = .exe
else
  ifeq ($(systype),Linux)  
    ifeq ($(x86compiler),gcc) 
      #Compiler information
      CC = gcc
      OPTFLAGS = -O6
      COPTIONS = -DLINUX -D_FILE_OFFSET_BITS=64 -Wall -std=c99 -pedantic

      #Linker information
      LDOPTIONS =
      LD = gcc
    else
      #Compiler information
      CC = icc
      OPTFLAGS = -O3
      COPTIONS = -DLINUX -D_FILE_OFFSET_BITS=64 -std=c99 

      #Linker information
      LDOPTIONS =
      LD = icc

      ifeq ($(openmp),yes)
        COPTIONS  += -D__OPENMP__ -openmp -openmp-report2
        LDOPTIONS += -openmp
      endif
    endif
  endif


  ifeq ($(systype),Darwin)  
    #Compiler information
    CC = gcc
    OPTFLAGS = -O6
    COPTIONS = -DDARWIN -D_FILE_OFFSET_BITS=64 -Wall -std=c99 -pedantic

    #Linker information
    LDOPTIONS =
    LD = gcc

    #Linker information
    LDOPTIONS =
    LD = gcc
  endif


  ifeq ($(systype),SunOS)
    #Compiler information
    CC = /opt/SUNWspro/bin/cc
    OPTFLAGS = -xO4
    COPTIONS = -DSUNOS 
  
    #Linker information
    LDOPTIONS = 
    LD = /opt/SUNWspro/bin/cc
  endif

  ifeq ($(systype),CYGWIN)
    #Compiler information
    CC = gcc 
    OPTFLAGS = -O6
    COPTIONS = -DCYGWIN -DWIN32 -D_FILE_OFFSET_BITS=64 -Wall -std=c99 -pedantic -mno-cygwin 

    #Linker information
    LDOPTIONS = -mno-cygwin 
    LD = gcc
  endif


  #-------------------------------------------------------------------
  # These defs are common among the GNU/GCC based systems
  #-------------------------------------------------------------------
  #Library creation information
  AR = ar crv $@ 
  RANLIB = ar -ts $@

  #Compile input/output file specification
  SOURCEFILE = -c $<
  OUTPUTFILE = -o $@

  # Standard file extensions
  OBJEXT = .o
  LIBEXT = .a
  EXEEXT = 
endif


#**************************************************************************
DMALLOCINC = 
DMALLOCFLAGS =
DEBUGFLAGS =
XTRAFLAGS = 

ifeq ($(dmalloc),yes)
  DMALLOCINC = -I$(HOME)/local/include
  DMALLOCFLAGS = -DDMALLOC
  OPTFLAGS = -g
endif

ifeq ($(debug),yes)
  DEBUGFLAGS = -DDEBUG
  OPTFLAGS = -g
endif

ifeq ($(gdb),yes)
  OPTFLAGS += -g
endif

ifeq ($(pcre),yes)
  XTRAFLAGS += -D__WITHPCRE__
endif
#**************************************************************************



#**************************************************************************
INCLUDES = -I./ $(INCDIR) $(DMALLOCINC)
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(XTRAFLAGS) $(DMALLOCFLAGS) $(DEBUGFLAGS) $(INCLUDES)
#**************************************************************************


#**************************************************************************
# Create the build directory if it does not exist
#**************************************************************************
BUILDDIR = ../builds/$(systype)-$(cputype)
$(shell mkdir -p $(BUILDDIR))

OBJS    = $(patsubst %.c, $(BUILDDIR)/%$(OBJEXT), $(wildcard *.c))
HEADERS = $(wildcard *.h) Makefile


$(BUILDDIR)/libGKlib$(LIBEXT): $(OBJS)
	$(AR) $(OBJS)
	$(RANLIB) 

doc:
	doxygen ../doxygen/Doxyfile

pubdoc:
	doxygen ../doxygen/Doxyfile
	scp -r ../doxygen/html/* dminers.dtc.umn.edu:/var/www/wwwroot/dminers/files/doxy/GKlib 
	ssh dminers.dtc.umn.edu "chmod -R a+rX /var/www/wwwroot/dminers/files/doxy/GKlib"

clean:
	rm -f $(OBJS)

realclean:
	rm -f $(OBJS) ; rm -f $(BUILDDIR)/libGKlib$(LIBEXT)


$(OBJS) : $(HEADERS) Makefile


$(BUILDDIR)/%$(OBJEXT) : %.c
	$(CC) $(CFLAGS) $(SOURCEFILE) $(OUTPUTFILE)


