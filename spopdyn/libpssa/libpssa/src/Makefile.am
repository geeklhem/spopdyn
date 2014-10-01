
LINKLIBS = $(OPT_LIBS) $(GSL_LIBS)

lib_LTLIBRARIES = libpssa.la
libpssa_la_SOURCES = PSSA.cpp \
			datamodel/DataModels.cpp \
			datamodel/PSSACR_Bins.cpp \
			datamodel/SimulationInfo.cpp \
			grouping/GroupingModule.cpp \
			grouping/GroupingModule_DM.cpp \
			grouping/GroupingModule_PDM.cpp \
			grouping/GroupingModule_PSSACR.cpp \
			grouping/GroupingModule_SPDM.cpp \
			sampling/CompositionRejectionSampler.cpp \
			sampling/SamplingModule.cpp \
			sampling/SamplingModule_DM.cpp \
			sampling/SamplingModule_PDM.cpp \
			sampling/SamplingModule_PSSACR.cpp \
			sampling/SamplingModule_SPDM.cpp \
			sampling/sha1.cpp \
			sampling/sha1.h \
			update/UpdateModule.cpp \
			update/UpdateModule_DM.cpp \
			update/UpdateModule_PDM.cpp \
			update/UpdateModule_PSSACR.cpp \
			update/UpdateModule_SPDM.cpp \
			gnuplot/IBgnuplot.cpp \
			vtk/vtk.cpp \
			MPI/mpi_spread.cpp \
			../dependencies/dSFMT-src-2.1/dSFMT-params.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params11213.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params1279.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params132049.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params19937.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params216091.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params2203.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params4253.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params44497.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params521.h \
			../dependencies/dSFMT-src-2.1/dSFMT-params86243.h \
			../dependencies/dSFMT-src-2.1/dSFMT.c \
			../dependencies/dSFMT-src-2.1/dSFMT.h \
			../dependencies/dSFMT-src-2.1/gsl_dSFMT.h


include_HEADERS = ../include/datamodel/CompositionRejectionSamplerData.h \
			../include/datamodel/DataModel.h \
			../include/datamodel/DataModel_DM.h \
			../include/datamodel/DataModel_PDM.h \
			../include/datamodel/DataModel_PSSACR.h \
			../include/datamodel/DataModel_SPDM.h \
			../include/datamodel/PSSACR_Bins.h \
			../include/datamodel/SimulationInfo.h \
			../include/grouping/GroupingModule.h \
			../include/grouping/GroupingModule_DM.h \
			../include/grouping/GroupingModule_PDM.h \
			../include/grouping/GroupingModule_PSSACR.h \
			../include/grouping/GroupingModule_SPDM.h \
			../include/sampling/CompositionRejectionSampler.h \
			../include/sampling/SamplingModule.h \
			../include/sampling/SamplingModule_DM.h \
			../include/sampling/SamplingModule_PDM.h \
			../include/sampling/SamplingModule_PSSACR.h \
			../include/sampling/SamplingModule_SPDM.h \
			../include/gnuplot/IBgnuplot.h \
			../include/update/UpdateModule.h \
			../include/update/UpdateModule_DM.h \
			../include/update/UpdateModule_PDM.h \
			../include/update/UpdateModule_PSSACR.h \
			../include/update/UpdateModule_SPDM.h \
			../include/vtk/vtk.h \
			../include/MPI/mpi_spread.h \
			../include/PSSA.h \
			../include/optimization.h \
			../include/oscompat.h \
			../include/sampling \
			../include/stdheaders.h \
			../include/typedefs.h



libpssa_la_CFLAGS = -DUNIX -rdynamic $(GSL_CFLAGS)
libpssa_la_CXXFLAGS = -DUNIX -rdynamic $(GSL_CFLAGS)
libpssa_la_LIBADD = $(LINKLIBS)
