 

LINKLIBS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) -L../libpssa/src/.libs -lpssa
		
bin_PROGRAMS = pssa
pssa_SOURCES = main.cpp
pssa_CFLAGS = -DUNIX -rdynamic -I../libpssa/include/ $(BOOST_CPPFLAGS)
pssa_CXXFLAGS = -DUNIX -rdynamic -I../libpssa/include/ $(BOOST_CPPFLAGS)
pssa_LDADD = $(LINKLIBS)
