LINKLIBS = $(OPT_LIBS) -L../../libpssa/src/.libs -lpssa

bin_PROGRAMS = grays_ms grays_alpha grays_sp grays_alpha7 grays_sp7 grays_mv7
grays_ms_SOURCES = main_ms.cpp
grays_sp_SOURCES = main_sp.cpp
grays_alpha_SOURCES = main_alpha.cpp
grays_alpha7_SOURCES = main_alpha7.cpp
grays_sp7_SOURCES = main_sp7.cpp
grays_mv7_SOURCES = main_mv7.cpp

graysdir=.
nobase_dist_grays_DATA = sbml/grayscott_mv.sbml \
sbml/grayscott_sp.sbml \
sbml/grayscott_alpha.sbml \
sbml/grayscott_alpha7.sbml \
sbml/grayscott_sp7.sbml \
sbml/grayscott_mv7.sbml \
out/out.txt

grays_ms_CFLAGS = -DUNIX -rdynamic -I../../libpssa/include $(GSL_CFLAGS)
grays_ms_CXXFLAGS = -DUNIX -rdynamic  -I../../libpssa/include $(GSL_CFLAGS)
grays_ms_LDADD = $(LINKLIBS)

grays_sp_CFLAGS = -DUNIX -rdynamic -I../../libpssa/include $(GSL_CFLAGS)
grays_sp_CXXFLAGS = -DUNIX -rdynamic  -I../../libpssa/include $(GSL_CFLAGS)
grays_sp_LDADD = $(LINKLIBS)

grays_alpha_CFLAGS = -DUNIX -rdynamic -I../../libpssa/include $(GSL_CFLAGS)
grays_alpha_CXXFLAGS = -DUNIX -rdynamic  -I../../libpssa/include $(GSL_CFLAGS)
grays_alpha_LDADD = $(LINKLIBS)

grays_alpha7_CFLAGS = -DUNIX -rdynamic -I../../libpssa/include $(GSL_CFLAGS)
grays_alpha7_CXXFLAGS = -DUNIX -rdynamic  -I../../libpssa/include $(GSL_CFLAGS)
grays_alpha7_LDADD = $(LINKLIBS)

grays_sp7_CFLAGS = -DUNIX -rdynamic -I../../libpssa/include $(GSL_CFLAGS)
grays_sp7_CXXFLAGS = -DUNIX -rdynamic  -I../../libpssa/include $(GSL_CFLAGS)
grays_sp7_LDADD = $(LINKLIBS)

grays_mv7_CFLAGS = -DUNIX -rdynamic -I../../libpssa/include $(GSL_CFLAGS)
grays_mv7_CXXFLAGS = -DUNIX -rdynamic  -I../../libpssa/include $(GSL_CFLAGS)
grays_mv7_LDADD = $(LINKLIBS)
