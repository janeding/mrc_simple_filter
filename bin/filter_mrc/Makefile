SHELL = /bin/sh

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cpp .o .a
.cpp.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


INTERNAL_LIB_PATH = ../../lib


LIBS = \
-L$(INTERNAL_LIB_PATH)/mrc_simple -lmrc_simple \
-lfftw3 \
-lm


LIB_FILES = \
$(INTERNAL_LIB_PATH)/mrc_simple/libmrc_simple.a


INCLUDES = \
-I$(INTERNAL_LIB_PATH)/err \
-I$(INTERNAL_LIB_PATH)/alloc \
-I$(INTERNAL_LIB_PATH)/filter \
-I$(INTERNAL_LIB_PATH)/threshold \
-I$(INTERNAL_LIB_PATH)/mrc_simple


OBJECT_FILES = settings.o filter_mrc.o


OBJECT_SRC = settings.cpp filter_mrc.cpp

RCS_FILES = $(OBJECT_SRC)



#Default target: the filter_mrc binary
#                "filter_mrc" is the lattice polymer simulator
filter_mrc: $(OBJECT_FILES) $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o filter_mrc \
	$(OBJECT_FILES) \
	$(LIBS)


install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" filter_mrc
	cp -f filter_mrc $(INSTALL_PATH)
	#cp -f docs_filter_mrc.txt $(INSTALL_PATH)

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core filter_mrc $(COMPILER_TEMP_FILES) $(LINKER_TEMP_FILES)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/filter_mrc
	#rm -f $(INSTALL_PATH)/docs_filter_mrc.txt

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
