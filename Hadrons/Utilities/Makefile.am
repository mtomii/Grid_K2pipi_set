bin_PROGRAMS = HadronsXmlRun HadronsXmlValidate HadronsFermionEP64To32 HadronsContractor HadronsContractorBenchmark K2pipiContractor BubbleContractor FlexibleContractor

HadronsXmlRun_SOURCES = HadronsXmlRun.cc
HadronsXmlRun_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

HadronsXmlValidate_SOURCES = HadronsXmlValidate.cc
HadronsXmlValidate_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

HadronsFermionEP64To32_SOURCES  = EigenPackCast.cc
HadronsFermionEP64To32_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
HadronsFermionEP64To32_LDADD    = ../libHadrons.a ../../Grid/libGrid.a

HadronsContractor_SOURCES = Contractor.cc
HadronsContractor_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

FlexibleContractor_SOURCES = FlexibleContractor.cc
FlexibleContractor_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

K2pipiContractor_SOURCES = K2pipiContractor.cc
K2pipiContractor_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

BubbleContractor_SOURCES = BubbleContractor.cc
BubbleContractor_LDADD   = ../libHadrons.a ../../Grid/libGrid.a

HadronsContractorBenchmark_SOURCES = ContractorBenchmark.cc
HadronsContractorBenchmark_LDADD   = ../libHadrons.a ../../Grid/libGrid.a
