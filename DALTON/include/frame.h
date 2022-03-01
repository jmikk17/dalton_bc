! frame.h
!  POTNUC : nuclear repulsion energy
!  DIPNUC : nuclear dipole moment for origin DIPORG(1:3)
!  QPNUC  : nuclear quadrupole moment (not used)
      REAL*8  POTNUC, DIPNUC
      COMMON /FRAME/ POTNUC, DIPNUC(3)
!    & , QPNUC(6)
