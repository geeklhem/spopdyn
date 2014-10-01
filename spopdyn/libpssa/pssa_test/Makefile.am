LINKLIBS = $(OPT_LIBS) -L../libpssa/src/.libs -lpssa

bin_PROGRAMS = pssa_test
pssa_test_SOURCES = main.cpp \
			TestBase.cpp \
			TestBase.h \
			TestDiffusion.cpp \
			TestDiffusion.h \
			TestReaction.cpp \
			TestReaction.h \
			TestReactionDiffusion.cpp \
			TestReactionDiffusion.h

pssa_testdir=.
nobase_dist_pssa_test_DATA = sbml/Diffusion.sbml \
		sbml/Multimerization.sbml

outdir=.
nobase_dist_out_DATA = out/out.txt

pssa_test_CFLAGS = -DUNIX -rdynamic -I../libpssa/include $(GSL_CFLAGS)
pssa_test_CXXFLAGS = -DUNIX -rdynamic  -I../libpssa/include $(GSL_CFLAGS)
pssa_test_LDADD = $(LINKLIBS)

