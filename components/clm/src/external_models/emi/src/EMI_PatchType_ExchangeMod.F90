module EMI_PatchType_Exchange
  
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use clm_varctl                            , only : iulog
  use ExternalModelInterfaceDataMod         , only : emi_data_list, emi_data
  use ExternalModelIntefaceDataDimensionMod , only : emi_data_dimension_list_type
  !
  implicit none
  !
  !
  public :: EMI_Pack_PatchType_for_EM

contains
!-----------------------------------------------------------------------
  subroutine EMI_Pack_PatchType_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's patch type for EM
    !
    ! !USES:
    use ExternalModelConstants , only : L2E_PATCH_ACTIVE
    use ExternalModelConstants , only : L2E_PATCH_TYPE
    use ExternalModelConstants , only : L2E_PATCH_WT_COL
    use VegetationType         , only : veg_pp                
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: fp,p
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_PATCH_ACTIVE)
             do fp = 1, num_filter
                p = filter(fp)
                if (veg_pp%active(p)) cur_data%data_int_1d(p) = 1
             enddo
             cur_data%is_set = .true.

          case (L2E_PATCH_TYPE)
             do fp = 1, num_filter
                p = filter(fp)
                cur_data%data_int_1d(p) = veg_pp%itype(p)
             enddo
             cur_data%is_set = .true.

          case (L2E_PATCH_WT_COL)
             do fp = 1, num_filter
                p = filter(fp)
                cur_data%data_real_1d(p) = veg_pp%wtcol(p)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMI_Pack_PatchType_for_EM

end module EMI_PatchType_Exchange
