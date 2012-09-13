 PSISRC = $(HOME)/Desktop/psi3
  PSIOBJ = /usr/local/psi
 
  VPATH = $(PSIOBJ)/lib
  PSILIBS = -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1
  CXXFLAGS = -g -Wall
  CPPFLAGS = -I$(PSIOBJ)/include
 
  PROGRAM = myscf
  PROGSRC = myscf.cc
  BINOBJ = $(PROGSRC:%.cc=%.o)
 
  $(PROGRAM): $(BINOBJ) $(PSILIBS) /usr/local/lib/ATLAS/objdir/lib/liblapack.a /usr/local/lib/ATLAS/objdir/lib/libcblas.a /usr/local/lib/ATLAS/objdir/lib/libf77blas.a /usr/local/lib/ATLAS/objdir/lib/libatlas.a /usr/lib/gcc/x86_64-linux-gnu/4.6/libgfortran.a 
	$(CXX) $(CXXFLAGS) $^ -o $@
 
  install: $(PROGRAM)
	install -c $(PROGRAM) $(HOME)/bin
 
  clean:
	/bin/rm -f $(PROGRAM) *.o core *~


