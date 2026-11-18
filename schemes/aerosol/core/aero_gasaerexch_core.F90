! Portable gas-aerosol exchange (gas condensation onto aerosol particles)
! science shared by the aerosol models (MAM, CARMA):
!  - partitioning of non-volatile gas uptake among modes/bins, reversible
!    uptake toward an equilibrium gas mixing ratio, and the nh4/so4 limit
!  - reversible (or irreversible) SOA gas-aerosol partitioning
!
! Model-specific bookkeeping (constituent indices, renaming, aging, history
! output, ...) stays with the calling aerosol model. The H2SO4 gas-to-aerosol
! mass transfer (uptake) rates depend on the aerosol size representation and
! are provided by the host's concrete aerosol_properties classes
! (gas_uptkrate procedure).
module aero_gasaerexch_core
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aero_gasaerexch_partition
  public :: aero_gasaerexch_reversible
  public :: aero_gasaerexch_nh4_limit
  public :: aero_gasaerexch_soa
  public :: aero_gasaerexch_soa_column
  public :: aero_soaexch

  ! gas diffusivity ratios relative to h2so4
  ! (differences in fuchs-sutugin and accom coef ignored)
  real(kind_phys), parameter, public :: h2so4_diffus_ratio = 1.0_kind_phys ! h2so4/h2so4
  real(kind_phys), parameter, public :: soa_diffus_ratio = 0.81_kind_phys  ! soa/h2so4
  real(kind_phys), parameter, public :: nh3_diffus_ratio = 2.08_kind_phys  ! nh3/h2so4

  ! methods for soa gas uptake
  integer, parameter, public :: method_soa_none = 0         ! no uptake
  integer, parameter, public :: method_soa_irreversible = 1 ! irreversible uptake done like h2so4 uptake
  integer, parameter, public :: method_soa_reversible = 2   ! reversible uptake using aero_soaexch

contains

  ! partitions the uptake of a non-volatile gas among modes/bins in a grid cell:
  !   fgain(n) = fraction of the total gas uptake going to mode/bin n
  !   sum_uprt = total gas-to-aerosol mass transfer rate (1/s)
  !   avg_uprt = average uptake rate over the time step, such that
  !              (gas mixing ratio)*avg_uprt is the total uptake (mol/mol/s)
  pure subroutine aero_gasaerexch_partition(nbins, do_bin, uptkrate, diffus_ratio, deltatxx, &
                                            fgain, sum_uprt, avg_uprt)
    integer, intent(in) :: nbins                     ! number of modes/bins
    logical, intent(in) :: do_bin(:)                 ! (nbins) mode/bin takes up the gas
    real(kind_phys), intent(in) :: uptkrate(:)       ! (nbins) h2so4 gas-to-aerosol mass transfer rates [s-1]
    real(kind_phys), intent(in) :: diffus_ratio      ! gas diffusivity ratio relative to h2so4 [1]
    real(kind_phys), intent(in) :: deltatxx          ! time step, slightly increased [s]
    real(kind_phys), intent(out) :: fgain(:)         ! (nbins) fraction of the gas uptake going to each mode/bin [1]
    real(kind_phys), intent(out) :: sum_uprt         ! total gas-to-aerosol mass transfer rate [s-1]
    real(kind_phys), intent(out) :: avg_uprt         ! average uptake rate over the time step [s-1]

    integer :: n

    sum_uprt = 0.0_kind_phys
    do n = 1, nbins
      if (do_bin(n)) then
        fgain(n) = uptkrate(n)*diffus_ratio
        sum_uprt = sum_uprt + fgain(n)
      else
        fgain(n) = 0.0_kind_phys
      end if
    end do

    if (sum_uprt > 0.0_kind_phys) then
      do n = 1, nbins
        fgain(n) = fgain(n)/sum_uprt
      end do
    end if
    ! at this point (sum_uprt <= 0.0) only when all the fgain are zero

    ! uptake amount (fraction of gas uptaken) over deltat
    avg_uprt = (1.0_kind_phys - exp(-deltatxx*sum_uprt))/deltatxx

  end subroutine aero_gasaerexch_partition

  ! computes aerosol tendencies of a gas (e.g., h2so4) with reversible uptake
  ! in a grid cell: the gas relaxes toward an equilibrium mixing ratio given
  ! for each mode/bin (e.g., sulfeq above the tropopause)
  pure subroutine aero_gasaerexch_reversible(nbins, do_bin, uptkrate, gas_equil, deltatxx, &
                                             sum_uprt, gas_bgn, aer_bgn, dqdt_aer, sum_dqdt_aer)
    integer, intent(in) :: nbins                     ! number of modes/bins
    logical, intent(in) :: do_bin(:)                 ! (nbins) mode/bin takes up the gas
    real(kind_phys), intent(in) :: uptkrate(:)       ! (nbins) gas-to-aerosol mass transfer rates [s-1]
    real(kind_phys), intent(in) :: gas_equil(:)      ! (nbins) equilibrium gas mixing ratio over each mode/bin [mol mol-1]
    real(kind_phys), intent(in) :: deltatxx          ! time step, slightly increased [s]
    real(kind_phys), intent(in) :: sum_uprt          ! total gas-to-aerosol mass transfer rate [s-1]
    real(kind_phys), intent(in) :: gas_bgn           ! initial gas mixing ratio [mol mol-1]
    real(kind_phys), intent(in) :: aer_bgn(:)        ! (nbins) initial aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(inout) :: dqdt_aer(:)    ! (nbins) aerosol tendencies [mol mol-1 s-1];
                                                     ! left unchanged where .not.do_bin
    real(kind_phys), intent(out) :: sum_dqdt_aer     ! total aerosol tendency [mol mol-1 s-1]

    integer :: n
    real(kind_phys) :: tmp_kxt, tmp_pxt
    real(kind_phys) :: tmp_aer_end
    real(kind_phys) :: tmp_gas_avg, tmp_gas_equ

    tmp_kxt = deltatxx*sum_uprt  ! sum over modes of uptake_rate*deltat
    tmp_pxt = 0.0_kind_phys
    do n = 1, nbins
      if (.not. do_bin(n)) cycle
      tmp_pxt = tmp_pxt + uptkrate(n)*gas_equil(n)
    end do
    tmp_pxt = max(0.0_kind_phys, tmp_pxt*deltatxx)  ! sum over modes of uptake_rate*gas_equil*deltat
    ! calc avg gas over deltat
    if (tmp_kxt >= 1.0e-5_kind_phys) then
      ! exponential decay towards equilibrium value solution
      tmp_gas_equ = tmp_pxt/tmp_kxt
      tmp_gas_avg = tmp_gas_equ + (gas_bgn - tmp_gas_equ)*(1.0_kind_phys - exp(-tmp_kxt))/tmp_kxt
    else
      ! first order approx for tmp_kxt small
      tmp_gas_avg = gas_bgn*(1.0_kind_phys - 0.5_kind_phys*tmp_kxt) + 0.5_kind_phys*tmp_pxt
    end if
    sum_dqdt_aer = 0.0_kind_phys
    do n = 1, nbins
      if (.not. do_bin(n)) cycle
      ! calc change to aerosol in mode n
      tmp_aer_end = aer_bgn(n) + deltatxx*uptkrate(n)*(tmp_gas_avg - gas_equil(n))
      tmp_aer_end = max(0.0_kind_phys, tmp_aer_end)
      dqdt_aer(n) = (tmp_aer_end - aer_bgn(n))/deltatxx
      sum_dqdt_aer = sum_dqdt_aer + dqdt_aer(n)
    end do

  end subroutine aero_gasaerexch_reversible

  ! computes nh4 aerosol tendencies from nh3 gas uptake in a grid cell, with
  ! the nh4/so4 molar ratio limited to 2
  pure subroutine aero_gasaerexch_nh4_limit(nbins, fgain_nh4, sum_dqdt_nh4, deltat, deltatxx, &
                                            qold_nh4, qold_so4, dqdt_so4, dqdt_nh4, sum_dqdt_nh4_b)
    integer, intent(in) :: nbins                     ! number of modes/bins
    real(kind_phys), intent(in) :: fgain_nh4(:)      ! (nbins) fraction of the nh3 uptake going to each mode/bin [1]
    real(kind_phys), intent(in) :: sum_dqdt_nh4      ! total nh3 uptake without the limit [mol mol-1 s-1]
    real(kind_phys), intent(in) :: deltat            ! time step [s]
    real(kind_phys), intent(in) :: deltatxx          ! time step, slightly increased [s]
    real(kind_phys), intent(in) :: qold_nh4(:)       ! (nbins) initial nh4 aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(in) :: qold_so4(:)       ! (nbins) initial so4 aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(in) :: dqdt_so4(:)       ! (nbins) so4 aerosol tendencies [mol mol-1 s-1]
    real(kind_phys), intent(out) :: dqdt_nh4(:)      ! (nbins) nh4 aerosol tendencies [mol mol-1 s-1]
    real(kind_phys), intent(out) :: sum_dqdt_nh4_b   ! total nh3 uptake with the limit [mol mol-1 s-1]

    integer :: n
    real(kind_phys) :: qmax_nh4, qnew_nh4, qnew_so4

    sum_dqdt_nh4_b = 0.0_kind_phys
    do n = 1, nbins
      dqdt_nh4(n) = fgain_nh4(n)*sum_dqdt_nh4
      qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
      qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
      qmax_nh4 = 2.0_kind_phys*qnew_so4
      if (qnew_nh4 > qmax_nh4) then
        dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
      end if
      sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
    end do

  end subroutine aero_gasaerexch_nh4_limit

  ! computes soa gas uptake tendencies in a single grid cell
  !
  !   method_soa_reversible   -- reversible uptake using aero_soaexch
  !   method_soa_irreversible -- irreversible uptake done like h2so4 uptake
  !   method_soa_none         -- no uptake
  subroutine aero_gasaerexch_soa(method_soa, deltat, temp, pres, rgas, &
                                 nbins, nsoa, npoa, do_bin, do_soa, do_soag, uptkrate, &
                                 mw_poa_host, mw_soa_host, &
                                 qold_soag, qold_soa, qold_poa, &
                                 dqdt_soa, sum_dqdt_soa, errmsg, errflg)
    integer, intent(in) :: method_soa                   ! soa uptake method
    real(kind_phys), intent(in) :: deltat               ! time step [s]
    real(kind_phys), intent(in) :: temp                 ! air temperature [K]
    real(kind_phys), intent(in) :: pres                 ! air pressure [Pa]
    real(kind_phys), intent(in) :: rgas                 ! universal gas constant [J K-1 mol-1]
    integer, intent(in) :: nbins                        ! number of modes/bins
    integer, intent(in) :: nsoa                         ! number of soa species
    integer, intent(in) :: npoa                         ! number of poa species
    logical, intent(in) :: do_bin(:)                    ! (nbins) mode/bin takes part in reversible soa exchange
    logical, intent(in) :: do_soa(:, :)                 ! (nbins,nsoa) soa species takes up gas in mode/bin
    logical, intent(in) :: do_soag(:)                   ! (nsoa) soa gas species is present
    real(kind_phys), intent(in) :: uptkrate(:)          ! (nbins) h2so4 gas-to-aerosol mass transfer rates [s-1]
    real(kind_phys), intent(in) :: mw_poa_host(:)       ! (npoa) molec wght of poa used in host code [kg kmol-1]
    real(kind_phys), intent(in) :: mw_soa_host(:)       ! (nsoa) molec wght of soa used in host code [kg kmol-1]
    real(kind_phys), intent(in) :: qold_soag(:)         ! (nsoa) soa gas mixing ratios [mol mol-1]
    real(kind_phys), intent(in) :: qold_soa(:, :)       ! (nbins,nsoa) soa aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(in) :: qold_poa(:, :)       ! (nbins,npoa) poa aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(out) :: dqdt_soa(:, :)      ! (nbins,nsoa) soa aerosol tendencies [mol mol-1 s-1]
    real(kind_phys), intent(out) :: sum_dqdt_soa(:)     ! (nsoa) total soa gas uptake [mol mol-1 s-1]
                                                        ! -- negative of gas tendency
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer, parameter :: niter_max = 1000

    integer :: jsoa, n, niter
    real(kind_phys) :: avg_uprt_soa(nsoa)
    real(kind_phys) :: deltatxx
    real(kind_phys) :: dqdt_soag(nsoa)
    real(kind_phys) :: fgain_soa(nbins, nsoa)
    real(kind_phys) :: sum_uprt_soa(nsoa)
    real(kind_phys) :: uptkrate_soa(nbins, nsoa)

    errmsg = ''
    errflg = 0

    ! use this for tendency calcs to avoid generating very small negative values
    deltatxx = deltat*(1.0_kind_phys + 1.0e-15_kind_phys)

    !   fgain_soa(n,jsoa) = fraction of total soa gas uptake going to mode/bin n
    sum_uprt_soa(:) = 0.0_kind_phys
    do n = 1, nbins
      do jsoa = 1, nsoa
        if (do_soa(n, jsoa)) then
          fgain_soa(n, jsoa) = uptkrate(n)*soa_diffus_ratio
          sum_uprt_soa(jsoa) = sum_uprt_soa(jsoa) + fgain_soa(n, jsoa)
        else
          fgain_soa(n, jsoa) = 0.0_kind_phys
        end if
        uptkrate_soa(n, jsoa) = fgain_soa(n, jsoa)
      end do
    end do

    do jsoa = 1, nsoa
      if (sum_uprt_soa(jsoa) > 0.0_kind_phys) then
        do n = 1, nbins
          fgain_soa(n, jsoa) = fgain_soa(n, jsoa)/sum_uprt_soa(jsoa)
        end do
      end if
    end do

    !   uptake amount (fraction of gas uptaken) over deltat
    do jsoa = 1, nsoa
      avg_uprt_soa(jsoa) = (1.0_kind_phys - exp(-deltatxx*sum_uprt_soa(jsoa)))/deltatxx
    end do

    !   sum_dqdt_soa = soa_a tendency from soa gas uptake (mol/mol/s)
    do jsoa = 1, nsoa
      if (do_soag(jsoa)) then
        sum_dqdt_soa(jsoa) = qold_soag(jsoa)*avg_uprt_soa(jsoa)
      else
        sum_dqdt_soa(jsoa) = 0.0_kind_phys
      end if
    end do

    select case (method_soa)
    case (method_soa_reversible)
      !   compute tendencies for soa gas and soa interstitial aerosol
      !   using soa parameterization
      call aero_soaexch(deltat, temp, pres, rgas, &
                        niter, niter_max, nbins, do_bin, npoa, nsoa, &
                        mw_poa_host, mw_soa_host, &
                        qold_soag, qold_soa, qold_poa, uptkrate_soa, &
                        dqdt_soag, dqdt_soa, errmsg, errflg)
      if (errflg /= 0) return
      sum_dqdt_soa(:) = -dqdt_soag(:)

    case (method_soa_irreversible)
      !   compute tendencies for soa interstitial aerosol
      !   due to simple gas uptake
      do jsoa = 1, nsoa
        do n = 1, nbins
          dqdt_soa(n, jsoa) = fgain_soa(n, jsoa)*sum_dqdt_soa(jsoa)
        end do
      end do

    case (method_soa_none)
      dqdt_soa(:, :) = 0.0_kind_phys
      sum_dqdt_soa(:) = 0.0_kind_phys

    case default
      dqdt_soa(:, :) = 0.0_kind_phys
      sum_dqdt_soa(:) = 0.0_kind_phys
      write (errmsg, '(a,i10)') 'aero_gasaerexch_soa: invalid method_soa =', method_soa
      errflg = 1
    end select

  end subroutine aero_gasaerexch_soa

  ! computes soa gas uptake tendencies for all grid cells of a chunk
  ! (levels top_lev to nlev), using aero_gasaerexch_soa in each cell
  !
  ! In each cell, the soa aerosol mixing ratio of mode/bin n and soa species
  ! jsoa is soa(i,k,n,jsoa) where do_soa(n,jsoa) is true, and zero otherwise.
  subroutine aero_gasaerexch_soa_column(method_soa, ncol, nlev, top_lev, deltat, temp, pmid, rgas, &
                                        nbins, nsoa, npoa, do_bin, do_soa, do_soag, uptkrate, &
                                        mw_poa_host, mw_soa_host, soag, soa, poa, &
                                        dqdt_soa, sum_dqdt_soa, errmsg, errflg)
    integer, intent(in) :: method_soa                   ! soa uptake method
    integer, intent(in) :: ncol                         ! number of columns
    integer, intent(in) :: nlev                         ! number of vertical levels
    integer, intent(in) :: top_lev                      ! top level index for computation
    real(kind_phys), intent(in) :: deltat               ! time step [s]
    real(kind_phys), intent(in) :: temp(:, :)           ! (ncol,nlev) air temperature [K]
    real(kind_phys), intent(in) :: pmid(:, :)           ! (ncol,nlev) air pressure [Pa]
    real(kind_phys), intent(in) :: rgas                 ! universal gas constant [J K-1 mol-1]
    integer, intent(in) :: nbins                        ! number of modes/bins
    integer, intent(in) :: nsoa                         ! number of soa species
    integer, intent(in) :: npoa                         ! number of poa species
    logical, intent(in) :: do_bin(:)                    ! (nbins) mode/bin takes part in reversible soa exchange
    logical, intent(in) :: do_soa(:, :)                 ! (nbins,nsoa) soa species takes up gas in mode/bin
    logical, intent(in) :: do_soag(:)                   ! (nsoa) soa gas species is present
    real(kind_phys), intent(in) :: uptkrate(:, :, :)    ! (ncol,nlev,nbins) h2so4 gas-to-aerosol mass transfer rates [s-1]
    real(kind_phys), intent(in) :: mw_poa_host(:)       ! (npoa) molec wght of poa used in host code [kg kmol-1]
    real(kind_phys), intent(in) :: mw_soa_host(:)       ! (nsoa) molec wght of soa used in host code [kg kmol-1]
    real(kind_phys), intent(in) :: soag(:, :, :)        ! (ncol,nlev,nsoa) soa gas mixing ratios [mol mol-1]
    real(kind_phys), intent(in) :: soa(:, :, :, :)      ! (ncol,nlev,nbins,nsoa) soa aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(in) :: poa(:, :, :, :)      ! (ncol,nlev,nbins,npoa) poa aerosol mixing ratios [mol mol-1]
    real(kind_phys), intent(out) :: dqdt_soa(:, :, :, :) ! (ncol,nlev,nbins,nsoa) soa aerosol tendencies [mol mol-1 s-1]
    real(kind_phys), intent(out) :: sum_dqdt_soa(:, :, :) ! (ncol,nlev,nsoa) total soa gas uptake [mol mol-1 s-1]
                                                          ! -- negative of gas tendency
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i, k, n, j
    real(kind_phys) :: uptkratebb(nbins)
    real(kind_phys) :: qold_soag(nsoa)
    real(kind_phys) :: qold_soa(nbins, nsoa)
    real(kind_phys) :: qold_poa(nbins, npoa)
    real(kind_phys) :: dqdt_soa_cell(nbins, nsoa)
    real(kind_phys) :: sum_dqdt_soa_cell(nsoa)

    errmsg = ''
    errflg = 0

    dqdt_soa(:, :, :, :) = 0.0_kind_phys
    sum_dqdt_soa(:, :, :) = 0.0_kind_phys

    do k = top_lev, nlev
      do i = 1, ncol
        do n = 1, nbins
          uptkratebb(n) = uptkrate(i, k, n)
          do j = 1, nsoa
            if (do_soa(n, j)) then
              qold_soa(n, j) = soa(i, k, n, j)
            else
              qold_soa(n, j) = 0.0_kind_phys
            end if
          end do
          do j = 1, npoa
            qold_poa(n, j) = poa(i, k, n, j)
          end do
        end do
        do j = 1, nsoa
          qold_soag(j) = soag(i, k, j)
        end do

        call aero_gasaerexch_soa(method_soa, deltat, temp(i, k), pmid(i, k), rgas, &
                                 nbins, nsoa, npoa, do_bin, do_soa, do_soag, uptkratebb, &
                                 mw_poa_host, mw_soa_host, &
                                 qold_soag, qold_soa, qold_poa, &
                                 dqdt_soa_cell, sum_dqdt_soa_cell, errmsg, errflg)
        if (errflg /= 0) return

        dqdt_soa(i, k, :, :) = dqdt_soa_cell(:, :)
        sum_dqdt_soa(i, k, :) = sum_dqdt_soa_cell(:)
      end do
    end do

  end subroutine aero_gasaerexch_soa_column

  ! calculates condensation/evaporation of "soa gas"
  ! to/from multiple aerosol modes or bins in 1 grid cell
  !
  ! key assumptions
  ! (1) ambient equilibrium vapor pressure of soa gas
  !     is given by p0_soa_298 and delh_vap_soa
  ! (2) equilibrium vapor pressure of soa gas at aerosol
  !     particle surface is given by raoults law in the form
  !     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
  ! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
  !
  ! Author: R. Easter and R. Zaveri
  ! Additions to run with multiple BC, SOA and POM's: Shrivastava et al., 2015
  ! 10.1002/2014JD022563
  subroutine aero_soaexch(dtfull, temp, pres, rgas, &
                          niter, niter_max, nbins, do_bin, ntot_poaspec, ntot_soaspec, &
                          mw_poa_host, mw_soa_host, &
                          g_soa_in, a_soa_in, a_poa_in, xferrate_in, &
                          g_soa_tend, a_soa_tend, errmsg, errflg)
    real(kind_phys), intent(in) :: dtfull               ! full integration time step [s]
    real(kind_phys), intent(in) :: temp                 ! air temperature [K]
    real(kind_phys), intent(in) :: pres                 ! air pressure [Pa]
    real(kind_phys), intent(in) :: rgas                 ! universal gas constant [J K-1 mol-1]
    integer, intent(out) :: niter                       ! number of iterations performed
    integer, intent(in) :: niter_max                    ! max allowed number of iterations
    integer, intent(in) :: nbins                        ! number of modes or bins
    logical, intent(in) :: do_bin(:)                    ! (nbins) mode/bin takes part in soa exchange
    integer, intent(in) :: ntot_poaspec                 ! number of poa species
    integer, intent(in) :: ntot_soaspec                 ! number of soa species
    real(kind_phys), intent(in) :: mw_poa_host(:)       ! (ntot_poaspec) molec wght of poa used in host code
    real(kind_phys), intent(in) :: mw_soa_host(:)       ! (ntot_soaspec) molec wght of soa used in host code
    real(kind_phys), intent(in) :: g_soa_in(:)          ! (ntot_soaspec) initial soa gas mixrat (mol/mol at host mw)
    real(kind_phys), intent(in) :: a_soa_in(:, :)       ! (nbins,ntot_soaspec) initial soa aerosol mixrat (mol/mol at host mw)
    real(kind_phys), intent(in) :: a_poa_in(:, :)       ! (nbins,ntot_poaspec) initial poa aerosol mixrat (mol/mol at host mw)
    real(kind_phys), intent(in) :: xferrate_in(:, :)    ! (nbins,ntot_soaspec) gas-aerosol mass transfer rate (1/s)
    real(kind_phys), intent(out) :: g_soa_tend(:)       ! (ntot_soaspec) soa gas mixrat tendency (mol/mol/s at host mw)
    real(kind_phys), intent(out) :: a_soa_tend(:, :)    ! (nbins,ntot_soaspec) soa aerosol mixrat tendency
                                                        ! (mol/mol/s at host mw)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: ll
    integer :: m

    real(kind_phys), parameter :: a_min1 = 1.0e-20_kind_phys
    real(kind_phys), parameter :: g_min1 = 1.0e-20_kind_phys
    real(kind_phys), parameter :: alpha = 0.05_kind_phys        ! parameter used in calc of time step
    real(kind_phys), parameter :: dtsub_fixed = -1.0_kind_phys  ! fixed sub-step for time integration (s)

    real(kind_phys) :: a_ooa_sum_tmp(nbins)          ! total ooa (=soa+opoa) in a mode/bin
    real(kind_phys) :: a_opoa(nbins)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
    real(kind_phys) :: a_soa(nbins, ntot_soaspec)     ! soa aerosol mixrat (mol/mol at actual mw)
    real(kind_phys) :: a_soa_tmp(nbins, ntot_soaspec) ! temporary soa aerosol mixrat (mol/mol)
    real(kind_phys) :: beta(nbins, ntot_soaspec)      ! dtcur*xferrate
    real(kind_phys) :: delh_vap_soa(ntot_soaspec)    ! heat of vaporization for gas soa (J/mol)
    real(kind_phys) :: del_g_soa_tmp(ntot_soaspec)
    real(kind_phys) :: dtcur                         ! current time step (s)
    real(kind_phys) :: dtmax                         ! = (dtfull-tcur)
    real(kind_phys) :: g0_soa(ntot_soaspec)          ! ambient soa gas equilib mixrat (mol/mol at actual mw)
    real(kind_phys) :: g_soa(ntot_soaspec)           ! soa gas mixrat (mol/mol at actual mw)
    real(kind_phys) :: g_star(nbins, ntot_soaspec)    ! soa gas mixrat that is in equilib
                                                     ! with each aerosol mode/bin (mol/mol)
    real(kind_phys) :: mw_poa(ntot_poaspec)          ! actual molec wght of poa
    real(kind_phys) :: mw_soa(ntot_soaspec)          ! actual molec wght of soa
    real(kind_phys) :: opoa_frac(ntot_poaspec)       ! fraction of poa that is opoa
    real(kind_phys) :: phi(nbins, ntot_soaspec)       ! "relative driving force"
    real(kind_phys) :: p0_soa(ntot_soaspec)          ! soa gas equilib vapor presssure (atm)
    real(kind_phys) :: p0_soa_298(ntot_soaspec)      ! soa gas equilib vapor presssure (atm) at 298 k
    real(kind_phys) :: sat(nbins, ntot_soaspec)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                     !    used by the numerical integration scheme -- it is not a saturation rato!
    real(kind_phys) :: tcur                          ! current integration time (from 0 s)
    real(kind_phys) :: tmpa, tmpb, tmpf
    real(kind_phys) :: tot_soa(ntot_soaspec)         ! g_soa + sum( a_soa(:) )
    real(kind_phys) :: xferrate(nbins, ntot_soaspec)  ! gas-aerosol mass transfer rate (1/s)

    errmsg = ''
    errflg = 0

    g_soa_tend(:) = 0.0_kind_phys
    a_soa_tend(:, :) = 0.0_kind_phys

    ! POA does not form solution with SOA for all runs; set opoa_frac=0.0_kind_phys  by Manish Shrivastava
    opoa_frac(:) = 0.0_kind_phys
    mw_poa(:) = 250.0_kind_phys
    mw_soa(:) = 250.0_kind_phys

    call soa_vbs_params(ntot_soaspec, p0_soa_298, delh_vap_soa, errmsg, errflg)
    if (errflg /= 0) return

    xferrate(:, :) = 0.0_kind_phys

    ! modes/bins involved in the soa gas-aerosol transfer
    ! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
    do m = 1, nbins
      if (.not. do_bin(m)) cycle
      do ll = 1, ntot_soaspec
        xferrate(m, ll) = xferrate_in(m, ll)
      end do
    end do

    ! convert incoming mixing ratios from mol/mol at the "host-code" molec. weight
    !    to mol/mol at the "actual" molec. weight (currently assumed to be 250.0)
    ! also
    !    force things to be non-negative
    !    calc tot_soa(ll)
    do ll = 1, ntot_soaspec
      tmpf = mw_soa_host(ll)/mw_soa(ll)
      g_soa(ll) = max(g_soa_in(ll), 0.0_kind_phys)*tmpf
      tot_soa(ll) = g_soa(ll)
      do m = 1, nbins
        if (.not. do_bin(m)) cycle
        a_soa(m, ll) = max(a_soa_in(m, ll), 0.0_kind_phys)*tmpf
        tot_soa(ll) = tot_soa(ll) + a_soa(m, ll)
      end do
    end do

    ! calc a_opoa
    ! note: a_opoa is overwritten for each poa species (not summed), and with
    ! opoa_frac=0 it is always zero; kept as is to preserve answers
    do m = 1, nbins
      if (.not. do_bin(m)) cycle
      a_opoa(m) = 0.0_kind_phys
      do ll = 1, ntot_poaspec
        a_opoa(m) = opoa_frac(ll)*a_poa_in(m, ll)
      end do
    end do

    ! calc ambient equilibrium soa gas
    do ll = 1, ntot_soaspec
      p0_soa(ll) = p0_soa_298(ll)* &
                   exp(-(delh_vap_soa(ll)/rgas)*((1.0_kind_phys/temp) - (1.0_kind_phys/298.0_kind_phys)))
      g0_soa(ll) = 1.01325e5_kind_phys*p0_soa(ll)/pres
    end do

    niter = 0
    tcur = 0.0_kind_phys
    dtcur = 0.0_kind_phys
    phi(:, :) = 0.0_kind_phys
    g_star(:, :) = 0.0_kind_phys

    ! integration loop -- does multiple substeps to reach dtfull
    time_loop: do while (tcur < dtfull - 1.0e-3_kind_phys)

      niter = niter + 1
      if (niter > niter_max) exit

      tmpa = 0.0_kind_phys  ! time integration parameter for all soa species
      do m = 1, nbins
        if (.not. do_bin(m)) cycle
        a_ooa_sum_tmp(m) = a_opoa(m) + sum(a_soa(m, 1:ntot_soaspec))
      end do
      do ll = 1, ntot_soaspec
        tmpb = 0.0_kind_phys  ! time integration parameter for a single soa species
        do m = 1, nbins
          if (.not. do_bin(m)) cycle
          sat(m, ll) = g0_soa(ll)/max(a_ooa_sum_tmp(m), a_min1)
          g_star(m, ll) = sat(m, ll)*a_soa(m, ll)
          phi(m, ll) = (g_soa(ll) - g_star(m, ll))/max(g_soa(ll), g_star(m, ll), g_min1)
          tmpb = tmpb + xferrate(m, ll)*abs(phi(m, ll))
        end do
        tmpa = max(tmpa, tmpb)
      end do

      if (dtsub_fixed > 0.0_kind_phys) then
        dtcur = dtsub_fixed
        tcur = tcur + dtcur
      else
        dtmax = dtfull - tcur
        if (dtmax*tmpa <= alpha) then
          ! here alpha/tmpa >= dtmax, so this is final substep
          dtcur = dtmax
          tcur = dtfull
        else
          dtcur = alpha/tmpa
          tcur = tcur + dtcur
        end if
      end if

      ! step 1 - for modes where soa is condensing, estimate "new" a_soa(m,ll)
      !    using an explicit calculation with "old" g_soa
      !    and g_star(m,ll) calculated using "old" a_soa(m,ll)
      ! do this to get better estimate of "new" a_soa(m,ll) and sat(m,ll)
      do m = 1, nbins
        if (.not. do_bin(m)) cycle
        do ll = 1, ntot_soaspec
          ! first ll loop calcs a_soa_tmp(m,ll) & a_ooa_sum_tmp
          a_soa_tmp(m, ll) = a_soa(m, ll)
          beta(m, ll) = dtcur*xferrate(m, ll)
          del_g_soa_tmp(ll) = g_soa(ll) - g_star(m, ll)
          if (del_g_soa_tmp(ll) > 0.0_kind_phys) then
            a_soa_tmp(m, ll) = a_soa(m, ll) + beta(m, ll)*del_g_soa_tmp(ll)
          end if
        end do
        a_ooa_sum_tmp(m) = a_opoa(m) + sum(a_soa_tmp(m, 1:ntot_soaspec))
        do ll = 1, ntot_soaspec
          ! second ll loop calcs sat & g_star
          if (del_g_soa_tmp(ll) > 0.0_kind_phys) then
            sat(m, ll) = g0_soa(ll)/max(a_ooa_sum_tmp(m), a_min1)
            g_star(m, ll) = sat(m, ll)*a_soa_tmp(m, ll)   ! this just needed for diagnostics
          end if
        end do
      end do

      ! step 2 - implicit in g_soa and semi-implicit in a_soa,
      !    with g_star(m,ll) calculated semi-implicitly
      do ll = 1, ntot_soaspec
        tmpa = 0.0_kind_phys
        tmpb = 0.0_kind_phys
        do m = 1, nbins
          if (.not. do_bin(m)) cycle
          tmpa = tmpa + a_soa(m, ll)/(1.0_kind_phys + beta(m, ll)*sat(m, ll))
          tmpb = tmpb + beta(m, ll)/(1.0_kind_phys + beta(m, ll)*sat(m, ll))
        end do

        g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_kind_phys + tmpb)
        g_soa(ll) = max(0.0_kind_phys, g_soa(ll))
        do m = 1, nbins
          if (.not. do_bin(m)) cycle
          a_soa(m, ll) = (a_soa(m, ll) + beta(m, ll)*g_soa(ll))/ &
                         (1.0_kind_phys + beta(m, ll)*sat(m, ll))
        end do
      end do

    end do time_loop

    ! calculate outgoing tendencies (at the host-code molec. weight)
    ! (a_soa & g_soa are at actual mw, but a_soa_in & g_soa_in are at host-code mw)
    do ll = 1, ntot_soaspec
      tmpf = mw_soa(ll)/mw_soa_host(ll)
      g_soa_tend(ll) = (g_soa(ll)*tmpf - g_soa_in(ll))/dtfull
      do m = 1, nbins
        if (.not. do_bin(m)) cycle
        a_soa_tend(m, ll) = (a_soa(m, ll)*tmpf - a_soa_in(m, ll))/dtfull
      end do
    end do

  end subroutine aero_soaexch

  ! volatility parameters of the soa species
  ! (New SOA properties added by Manish Shrivastava on 09/27/2012)
  pure subroutine soa_vbs_params(nsoa, p0_soa_298, delh_vap_soa, errmsg, errflg)
    integer, intent(in) :: nsoa                          ! number of soa species
    real(kind_phys), intent(out) :: p0_soa_298(:)        ! (nsoa) soa gas equilib vapor presssure at 298 k [atm]
    real(kind_phys), intent(out) :: delh_vap_soa(:)      ! (nsoa) heat of vaporization for gas soa [J mol-1]
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: j

    errmsg = ''
    errflg = 0

    select case (nsoa)
    case (1)
      p0_soa_298(:) = 9.7831E-11_kind_phys
      delh_vap_soa(:) = 131.0e3_kind_phys
    case (2)
      ! same for anthropogenic and biomass burning species
      p0_soa_298(:) = 1.0e-10_kind_phys
      delh_vap_soa(:) = 156.0e3_kind_phys
    case (5, 15)
      ! 5 volatility bins for each SOA class:
      !   nsoa=5  -- one combined class (biomass burning, fossil fuel, biogenic)
      !   nsoa=15 -- species 1-5 fossil fuel (anthropogenic), 6-10 biomass burning,
      !              11-15 biogenic SOA, based on Cappa et al., Reference needs to be updated
      ! delh_vap_soa following the numbers by Epstein et al., 2012
      ! For MW=250.0
      do j = 0, nsoa - 5, 5
        p0_soa_298(j + 1) = 9.7831E-13_kind_phys ! C*=0.01ug/m3
        p0_soa_298(j + 2) = 9.7831E-12_kind_phys ! C*=0.10ug/m3
        p0_soa_298(j + 3) = 9.7831E-11_kind_phys ! C*=1.0ug/m3
        p0_soa_298(j + 4) = 9.7831E-10_kind_phys ! C*=10.0ug/m3
        p0_soa_298(j + 5) = 9.7831E-9_kind_phys  ! C*=100.0ug/m3

        delh_vap_soa(j + 1) = 153.0e3_kind_phys
        delh_vap_soa(j + 2) = 142.0e3_kind_phys
        delh_vap_soa(j + 3) = 131.0e3_kind_phys
        delh_vap_soa(j + 4) = 120.0e3_kind_phys
        delh_vap_soa(j + 5) = 109.0e3_kind_phys
      end do
    case default
      p0_soa_298(:) = 0.0_kind_phys
      delh_vap_soa(:) = 0.0_kind_phys
      write (errmsg, '(a,i6)') 'aero_gasaerexch_core::soa_vbs_params: unsupported number of SOA species ', nsoa
      errflg = 1
    end select

  end subroutine soa_vbs_params

end module aero_gasaerexch_core
