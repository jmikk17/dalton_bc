! FILE: exeinf.h
!
!     ITRLVL_LAST : last integral transformation level (updated in sirtra.F)
!     LVLDRC_LAST : last Dirac integral transformation level (updated in sirtra.F)
!     *TRONV : variables used to control abatro.F(TROINV)
!     FTRCTL : true - force new integral transformation because AOs have changed (typically new geometry)
!     NEWCMO : true - new integral transformation needed because CMO has changed
!     (FT* true is abbreviation for FIRST call in a series)
!
      INTEGER         ITRLVL_LAST, LVLDRC_LAST
      LOGICAL         FTRONV, GTRONV, HTRONV, RTRONV,
     &                FTRCTL, NEWCMO
      COMMON /EXEINF/ ITRLVL_LAST, LVLDRC_LAST,
     &                FTRONV, GTRONV, HTRONV, RTRONV,
     &                FTRCTL, NEWCMO
! end of exeinf.h
