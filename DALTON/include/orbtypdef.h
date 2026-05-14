! orbtypdef.h
! IDBTYP positive values: 1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
!        negative values: involve frozen orbitals
      INTEGER     IDBTYP(4,4)
      CHARACTER*9 COBTYP(4)
      SAVE      IDBTYP, COBTYP
      DATA      IDBTYP/-1,-2,-3,-4,                                     &
     &                 -2, 1, 2, 4,                                     &
     &                 -3, 2, 3, 5,                                     &
     &                 -4, 4, 5, 6/
      DATA      COBTYP/'frozen   ','inactive ','active   ','secondary'/
