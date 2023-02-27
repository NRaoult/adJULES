!
! This source is/has been preprocessed - be sure to edit init_parms.F90
!



  SUBROUTINE reinit_parms







  USE ancil_info, ONLY : frac,halo_i,halo_j,ice_fract,     &
 &       ice_fract_ncat,land_index,land_pts,nice,ntiles,   &
 &       off_x,off_y,rows,row_length,tile_index,tile_pts,  &
 &       ssi_pts,sea_pts,sice_pts,sice_pts_ncat,           &
 &       ssi_index,sea_index,sice_index,sice_index_ncat,   &
 &       fssi,sea_frac,sice_frac,sice_frac_ncat

  USE switches, ONLY : can_model,routeOnly,ltimer,l_aggregate
  USE prognostics, ONLY : canht_ft,di,di_ncat,lai,tstar_tile
  USE p_s_parms, ONLY : catch,catch_snow,infil_tile,               &
 &            satcon,z0_tile                      
  USE fluxes, ONLY : tstar

  USE u_v_grid, ONLY : dtrdz_charney_grid_1
  USE coastal, ONLY : flandg,tstar_land,tstar_sea,tstar_sice,tstar_ssi
 
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: I,J,L,N         ! WORK Loop counters

!-----------------------------------------------------------------------
! If only doing routing, nothing in here need be set.
!-----------------------------------------------------------------------
  IF ( routeOnly ) RETURN

!-----------------------------------------------------------------------
! Set up ancillary information variables
!-----------------------------------------------------------------------

! Set halo information to zero (no halos needed for stand alone model)
  OFF_X=0
  OFF_Y=0
  HALO_I=0
  HALO_J=0

! Set logical ltimer to false as no timing information available
  LTIMER=.FALSE.

!-----------------------------------------------------------------------
! Calculate surface parameters.
!-----------------------------------------------------------------------
  CALL SPARM (LAND_PTS,NTILES,CAN_MODEL,L_AGGREGATE                     &
 &,                 TILE_PTS,TILE_INDEX,FRAC,CANHT_FT,LAI,SATCON        &
 &,                 CATCH_SNOW,CATCH,INFIL_TILE,Z0_TILE)

!-----------------------------------------------------------------------
! Set up index for sea and sea-ice
!-----------------------------------------------------------------------
  SSI_PTS = 0
  SSI_INDEX(:)=0
  DO J=1,ROWS
    DO I=1,ROW_LENGTH
      SSI_PTS=SSI_PTS + 1
      IF ( FLANDG(I,J) < 1.0 ) THEN
        SSI_INDEX(SSI_PTS) = (J-1)*ROW_LENGTH + I
      ENDIF
      FSSI(I,J)=1.0 - FLANDG(I,J)
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set sea ice fraction.
!-----------------------------------------------------------------------
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      ICE_FRACT(I,J)=0.0
      DI(I,J)=0.0
      DO N=1,NICE
        ICE_FRACT(I,J)=ICE_FRACT(I,J)+ICE_FRACT_NCAT(I,J,N)
        DI(I,J)=DI(I,J)+ICE_FRACT_NCAT(I,J,N)*DI_NCAT(I,J,N)
      ENDDO
    ENDDO
  ENDDO








!-----------------------------------------------------------------------
! Initialise sea and sea-ice indices.
!-----------------------------------------------------------------------
  SEA_PTS = 0
  SICE_PTS = 0
  SEA_INDEX(:)=0
  SICE_INDEX(:)=0
  SICE_FRAC(:)=0.0
  SEA_FRAC(:)=0.0
  DO L=1,SSI_PTS
    J=(SSI_INDEX(L)-1)/ROW_LENGTH + 1
    I = SSI_INDEX(L) - (J-1)*ROW_LENGTH
    IF ( SSI_INDEX(L) > 0 ) THEN
      IF ( ICE_FRACT(I,J) > 0.0 ) THEN
        SICE_PTS=SICE_PTS+1
        SICE_INDEX(SICE_PTS)=L
        SICE_FRAC(L)=ICE_FRACT(I,J)
      ENDIF
      IF( ICE_FRACT(I,J) < 1.0 ) THEN
        SEA_PTS=SEA_PTS+1
        SEA_INDEX(SEA_PTS)=L
        SEA_FRAC(L)=1.0 - SICE_FRAC(L)
      ENDIF
    ENDIF
  ENDDO

  SICE_PTS_NCAT(:)=0
  SICE_INDEX_NCAT(:,:)=0
  SICE_FRAC_NCAT(:,:)=0.0
  DO N=1,NICE
    DO L=1,SSI_PTS
      J=(SSI_INDEX(L)-1)/ROW_LENGTH + 1
      I = SSI_INDEX(L) - (J-1)*ROW_LENGTH
      IF ( SSI_INDEX(L) > 0 ) THEN
        IF ( ICE_FRACT_NCAT(I,J,N) > 0.0 ) THEN
          SICE_PTS_NCAT(N)=SICE_PTS_NCAT(N)+1
          SICE_INDEX_NCAT(SICE_PTS_NCAT(N),N)=L
          SICE_FRAC_NCAT(L,N)=ICE_FRACT_NCAT(I,J,N)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set up gridbox "prognostics".
!-----------------------------------------------------------------------

  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      TSTAR(I,J)=0.0
      TSTAR_LAND(I,J)=0.0
      TSTAR_SSI(I,J)=0.0
    ENDDO
  ENDDO

  DO L=1,LAND_PTS
    J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
    I= LAND_INDEX(L) - (J-1)*ROW_LENGTH
    IF ( L_AGGREGATE ) THEN
      TSTAR_LAND(I,J) = TSTAR_TILE(L,1)
    ELSE
      DO N=1,NTILES
        TSTAR_LAND(I,J)=TSTAR_LAND(I,J) + FRAC(L,N)*TSTAR_TILE(L,N)
      ENDDO
    ENDIF
  ENDDO
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      TSTAR_SSI(I,J)=(1.0-ICE_FRACT(I,J))*TSTAR_SEA(I,J)  &
                   + ICE_FRACT(I,J)*TSTAR_SICE(I,J)
      TSTAR(I,J)=FLANDG(I,J)*TSTAR_LAND(I,J)  &
                   + (1.0-FLANDG(I,J))*TSTAR_SSI(I,J)
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set up information on U, V and T grids (assume that att grids are the same)
!-----------------------------------------------------------------------
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      DTRDZ_CHARNEY_GRID_1(I,J)=0.0
    ENDDO
  ENDDO

  RETURN



  END SUBROUTINE reinit_parms

  

