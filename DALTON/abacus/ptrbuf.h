! FILE: ptrbuf.h -- only used in DALTON/abacus/abaptr.F
!                -- depends on: maxorb.h
      INTEGER, parameter :: MAXCHN = 2000

      INTEGER         LASTAD, IADR
      INTEGER         MX1BUF, L1BUF,  MEMS,   MEMT,  MXABUF,
     &                LABUF,  LDAMAX, MX2BUF, L2BUF, NCHAIN, NBLOCK
      COMMON /PTRBUF/
     &                LASTAD(MAXCHN), IADR(MAXCHN),
     &                MX1BUF, L1BUF,  MEMS,   MEMT,  MXABUF,
     &                LABUF,  LDAMAX, MX2BUF, L2BUF, NCHAIN, NBLOCK
! end of ptrbuf.h
