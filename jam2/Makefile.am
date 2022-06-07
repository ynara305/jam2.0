SUBDIRS = hadrons xsection initcond collision interaction meanfield fluid

# top-level library files
noinst_LTLIBRARIES = libJAM.la
libJAM_la_SOURCES = JamAna.cxx JAM.cxx JamAna.h JAM.h
jam2includedir = $(includedir)/jam2
jam2include_HEADERS = JamAna.h JAM.h

# combine libraries
lib_LTLIBRARIES = libjam2.la
libjam2_la_SOURCES =
libjam2_la_LIBADD = \
 libJAM.la \
 hadrons/libhadrons.la \
 xsection/libxsection.la \
 initcond/libinitcond.la \
 interaction/libinteraction.la \
 collision/libcollision.la \
 fluid/libfluid.la \
 meanfield/libmeanfield.la
# Note on "-version-info CURRENT:REVISION:AGE":
#   According to the libtool versioning, every time when the interface is
#   changed (i.e., when there are essential changes of public interfaces in
#   *.h), we need to increment CURRENT and AGE and reset REVISION to 0.  We
#   need to reset AGE to 0 when there is any breaking changes in the interface.
#   When *.cxx is changed without changing *.h, we need to increment REVISION.
libjam2_la_LDFLAGS = -version-info 0:0:0

bin_PROGRAMS = jam testjam
jam_MAIN = main.cxx
jam_LDFLAGS = $(LDFLAGS) $(LOADLIB_LDFLAGS)
#jam_SOURCES = $(jam_MAIN) myhist.cxx
jam_SOURCES = $(jam_MAIN)
jam_LDADD = $(libjam2_la_LIBADD) $(LOADLIB_LIB) $(PYTHIALIBS) -ldl

testjam_MAIN = maintest.cxx
testjam_LDFLAGS = $(LDFLAGS) $(LOADLIB_LDFLAGS)
testjam_SOURCES = $(testjam_MAIN)
testjam_LDADD = $(jam_LDADD)

EXTRA_DIST = Manual.txt input
