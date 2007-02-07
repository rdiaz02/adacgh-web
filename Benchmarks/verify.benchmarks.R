## Simply verify we don't crash

fHMM(1, 2)
fCBS(1, 2)
## fBioHMM is very fragile ...
fBioHMM(1, 2)
fGLAD(1, 2)
mpiSetup(120)
fHMM.mpi(1, 2)
fCBS.mpi(1, 2)
fBioHMM.mpi(1, 2)
fGLAD.mpi(1, 2)
fACE.mpi(1, 2)
fCGHseg.mpi(1, 2)
fPSW.mpi(1, 2)
fWave.mpi(1, 2)
