# From the original Fermiqcd Library
#
# -DLINUX if on Linux/Unix NO -DLINUX on cygwin
# -DUSE_DOUBLE_PRECISION for double precision (the most efficient option)
# -DSSE2 -O3 for pentium4 optimizations
# -DPARALLEL is you have mpi then compile with mpiCC and run with mpirun 
#
#

#CC     = g++ -I../Libraries 
CC     = mpic++ -I../Libraries 
# Basic compilation flags (-O2 for speed)
#CFLAGS = -DLINUX -O3
# Flags for double precision and sse (P4 optiomizations)
#CFLAGS = -DLINUX -DSSE2 -O3
# Flags for parallel (mpi) double precision and sse
#CFLAGS = -DPARALLEL -DLINUX -DUSE_DOUBLE_PRECISION -DSSE2 -O2
# -Wno-write-strings suppresses depricated "conversion from string constant to ‘char*’" warning
CFLAGS = -DPARALLEL -DLINUX -O2 -g -Wno-write-strings


# Heat bath SU(N) with Polyakov loop measurements
pureSUN:: pureSUN.cpp ploop3.cpp
	${CC} ${CFLAGS} pureSUN.cpp -o pureSUN

# Heat bath SU(N) with Wilson flow
pureSUNwf_rk:: pureSUNwf_rk.cpp ploop3.cpp wilsonflow_rk.cpp readmilcascii.cpp
	${CC} ${CFLAGS} pureSUNwf_rk.cpp -o pureSUNwf_rk

pureSUNwflow:: pureSUNwflow.cpp ploop3.cpp wilsonflow_rk.cpp readmilcascii.cpp
	${CC} ${CFLAGS} pureSUNwflow.cpp -o pureSUNwflow


# Heat bath SU(N) with Wilson flow
pureSUNwf:: pureSUNwf.cpp ploop3.cpp wilsonflow.cpp
	${CC} ${CFLAGS} pureSUNwf.cpp -o pureSUNwf

pureSUNwf_1:: pureSUNwf_1.cpp ploop3.cpp wilsonflow.cpp
	${CC} ${CFLAGS} pureSUNwf_1.cpp -o pureSUNwf_1


