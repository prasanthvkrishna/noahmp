module NoahmpMosaicSortTileCatMod

!!! This module part of NoahMP Mosaic/Subgrid Tiling Scheme
!!! Purpose: To sort and identify most dominant lulc/soiltype/hydro types in a grid
!!!          identify dominant that contributed to > 90% grid area and 
!!!          scale the dominant ones to 100%

! ------------------------ Code history -----------------------------------
! Original code : Prasanth Valayamkunnath (IISER Thiruvananthapuram)
! Date          : July 10, 2025
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType

#ifdef MPP_LAND
  use module_mpp_land, only:mpp_land_bcast_int1, my_id, IO_id, mpp_status, numprocs, calculate_ntilemax_mpp
#endif

  implicit none

contains

!=== sort landuse/soiltype/hydrotype index based on area fraction and identify most dominant types

  subroutine NoahmpMosaicSortTileCat (NoahmpIO)

    implicit none

    type(NoahmpIO_type), intent(inout)  :: NoahmpIO

!---------------------------------------------------------------------
!  Local Variables
!---------------------------------------------------------------------

    integer                                           :: i, j, k, m 
    integer                                           :: max_index 
    integer                                           :: temp_idx
    integer                                           :: n_dominant
    integer,                allocatable, dimension(:) :: sorted_indices
    real(kind=kind_noahmp)                            :: cumulative_sum
    real(kind=kind_noahmp)                            :: temp_val
    real(kind=kind_noahmp), allocatable, dimension(:) :: frac_vec
    real(kind=kind_noahmp), allocatable, dimension(:) :: sorted_values
    real(kind=kind_noahmp), allocatable, dimension(:) :: rescaled_values
    integer                                           :: local_NTilesMax
    integer                                           :: global_NTilesMax
! -------------------------------------------------------------------------
    associate(                                                      &
              XSTART               =>  NoahmpIO%XSTART             ,&
              XEND                 =>  NoahmpIO%XEND               ,&
              YSTART               =>  NoahmpIO%YSTART             ,&
              YEND                 =>  NoahmpIO%YEND               ,&
              NumMosaicCat         =>  noahmpio%NumMosaicCat       ,&
              SubGrdFrac           =>  NoahmpIO%SubGrdFrac         ,&
              SubGrdFracRescaled   =>  NoahmpIO%SubGrdFracRescaled ,& 
              SubGrdIndexSorted    =>  NoahmpIO%SubGrdIndexSorted  ,&
              NumberOfTiles        =>  NoahmpIO%NumberOfTiles      ,&
              NTilesMax            =>  NoahmpIO%NTilesMax          ,&
              NTiles_user          =>  NoahmpIO%IOPT_MOSAIC_NTILES  &
             )
! -------------------------------------------------------------------------

!   Step1: Allocate local vectors
    if ( .not. allocated (sorted_indices) )     allocate ( sorted_indices     (NumMosaicCat) )
    if ( .not. allocated (frac_vec)       )     allocate ( frac_vec           (NumMosaicCat) )
    if ( .not. allocated (sorted_values)  )     allocate ( sorted_values      (NumMosaicCat) )   
    if ( .not. allocated (rescaled_values))     allocate ( rescaled_values    (NumMosaicCat) )


 ! initialize temp variables   
    temp_val         = 0.
    NTilesMax        = 1
    cumulative_sum   = 0.
    global_NTilesMax = 1
    local_NTilesMax  = 1

!   Step 2: Loop over each grid cell (i,j)
    do i = XSTART, XEND
       do j = YSTART, YEND

          ! Extract fractions and initialize sorting arrays
          do k = 1, NumMosaicCat
             if (k .eq. NoahmpIO%ISWATER) then     ! to skip water fractions in the grid
               SubGrdFrac(i,j,k) = 0.
             endif
             frac_vec(k)         = SubGrdFrac(i,j,k)
             sorted_indices(k)   = k
             sorted_values(k)    = frac_vec(k)
          end do

          ! Sort descending by fraction value (Selection Sort)
          do k = 1, NumMosaicCat - 1
            max_index = k
            do m = k + 1, NumMosaicCat
              if (sorted_values(m) > sorted_values(max_index)) then
                 max_index = m
              end if
            end do

            ! Swap values
            temp_val = sorted_values(k)
            sorted_values(k) = sorted_values(max_index)
            sorted_values(max_index) = temp_val

            temp_idx = sorted_indices(k)
            sorted_indices(k) = sorted_indices(max_index)
            sorted_indices(max_index) = temp_idx
          end do

          ! Step 3: Identify dominant types (sum > 0.9)
          cumulative_sum = 0.0
          n_dominant = 0
          do k = 1, NumMosaicCat
             if (sorted_values(k) == 0) exit
             ! Now limit the number of tiles to user defined number from the namelist
             if (k <= NTiles_user) then 
                n_dominant = k
                cumulative_sum = cumulative_sum + sorted_values(k)  !
                if (cumulative_sum > 0.9) exit
             else 
                exit
             endif
          end do
          ! Safety: If no valid tiles found, set at least 1 tile (water or default)
          if (n_dominant == 0) then
             n_dominant = 1
             sorted_indices(1) = NoahmpIO%ISWATER
             rescaled_values(1) = 1.0_kind_noahmp
          else
             ! Step 4: Rescale to sum = 1.0
             do k = 1, n_dominant
                rescaled_values(k) = sorted_values(k) / cumulative_sum
             end do
          endif

!          do k = 1, n_dominant
!             SubGrdFracRescaled(i,j,k)  =  rescaled_values(k)
!             SubGrdIndexSorted(i,j,k)   =  sorted_indices(k)
!             local_NTilesMax            =  min(max(NTilesMax,n_dominant), NTiles_user) ! maximum value across domain
!             NumberOfTiles(i,j)         =  min(local_NTilesMax,n_dominant) ! it is >0 and <= NTiles_user
!          enddo
!
!          if(i.eq.869.and.j.eq.21) then
!             print*,'NTilesMax=',NTilesMax
!             print*,'n_dominant=',n_dominant
!             print*,'NTiles_user=',NTiles_user
!             print*,'local_NTilesMax=',local_NTilesMax
!             print*,'NumberOfTiles(i,j)=',NumberOfTiles(i,j)
!          endif
!
          NumberOfTiles(i,j) = n_dominant 
          
          do k = 1, n_dominant
             SubGrdFracRescaled(i,j,k)  =  rescaled_values(k)
             SubGrdIndexSorted(i,j,k)   =  sorted_indices(k)
          enddo

          ! Track the running maximum across the ENTIRE local domain
          local_NTilesMax = max(local_NTilesMax, n_dominant)

      end do
    end do

!#ifdef MPP_LAND
!
!   call calculate_ntilemax_mpp(local_NTilesMax,global_NTilesMax)
!    ! Broadcast the final global max back to all ranks
!    call mpp_land_bcast_int1(global_NTilesMax)
!
!    ! Store final result
!    NTilesMax = min(global_NTilesMax, NTiles_user)
!#else
!    NTilesMax = local_NTilesMax
!#endif

#ifdef MPP_LAND
    call calculate_ntilemax_mpp(local_NTilesMax, global_NTilesMax)
    call mpp_land_bcast_int1(global_NTilesMax)
   
    NTilesMax = min(global_NTilesMax, NTiles_user)
#else
    NTilesMax = min(local_NTilesMax, NTiles_user)
#endif

    ! CRITICAL: Cap individual grid cell tile counts to the final allocated NTilesMax 
    ! to prevent downstream out-of-bounds errors if global max exceeded NTiles_user.
    do i = XSTART, XEND
      do j = YSTART, YEND
         NumberOfTiles(i,j) = min(NumberOfTiles(i,j), NTilesMax)
      enddo
    enddo
    print*,'NoahmpMosaicSortTileCat: ',"NTilesMax=",NTilesMax

    end associate

  end subroutine NoahmpMosaicSortTileCat

end module NoahmpMosaicSortTileCatMod
