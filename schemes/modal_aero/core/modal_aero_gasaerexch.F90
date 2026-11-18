! Portable code for modal aerosol gas-aerosol exchange.
! RCE 07.04.13:  Adapted from MIRAGE2 code
module modal_aero_gasaerexch
  use ccpp_kinds, only: kind_phys
  use aero_gasaerexch_core, only: aero_gasaerexch_partition, aero_gasaerexch_reversible, &
                                  aero_gasaerexch_nh4_limit, aero_gasaerexch_soa, &
                                  h2so4_diffus_ratio, nh3_diffus_ratio, method_soa_reversible

  implicit none
  private

  public :: modal_aero_gasaerexch_init
  public :: modal_aero_gasaerexch_run

  ! Primary-carbon aging (pcage) configuration: species are transferred from
  ! the primary-carbon mode (modefrm) to the accumulation mode (modetoo) when
  ! enough sulfate monolayers coat the particle.
  integer, protected, public :: maxspec_pcage     ! max number of species that can be aged

  integer, protected, public :: modefrm_pcage     ! source mode index for aging transfer
  integer, protected, public :: nspecfrm_pcage    ! number of species transferred during aging

  integer, protected, allocatable, public :: lspecfrm_pcage(:) ! pcnst indices of species in source mode
  integer, protected, allocatable, public :: lspectoo_pcage(:) ! pcnst indices of corresponding species in dest mode

  real(kind_phys), parameter, public :: n_so4_monolayers_pcage = 8.0_kind_phys

  ! number of so4(+nh4) monolayers needed to "age" a carbon particle
  ! thickness of the so4 monolayers (m)
  ! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
  !    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
  ! aging criterion is approximate so do not try to distinguish
  !    sulfuric acid, bisulfate, ammonium sulfate
  real(kind_phys), parameter, public :: &
    dr_so4_monolayers_pcage = n_so4_monolayers_pcage*4.76e-10_kind_phys

  ! this factor converts an soa volume to a volume of so4(+nh4)
  ! having same hygroscopicity as the soa
  real(kind_phys), protected, allocatable, public :: soa_equivso4_factor(:)

  ! Private module-level storage:

  ! Mode configuration
  integer :: ntot_amode_m, nsoa_m, npoa_m, nspec_max_m
  integer, allocatable :: nspec_amode_m(:)

  ! Species indices in pcnst-space (set by _init, converted to vmr-space in _run)
  integer :: idx_h2so4_m, idx_nh3_m, idx_msa_m
  integer, allocatable :: idx_soag_m(:)
  integer, allocatable :: idx_so4_a_m(:), idx_nh4_a_m(:)
  integer, allocatable :: idx_soa_a_m(:, :), idx_pom_a_m(:, :)
  integer, allocatable :: idx_num_m(:), idx_mass_m(:, :)

  ! Mode metadata
  real(kind_phys), allocatable :: alnsg_amode_m(:)
  real(kind_phys), allocatable :: specmw_amode_m(:, :), specdens_amode_m(:, :)

  ! Flags
  logical :: do_nh4g_m, do_msag_m, do_soag_any_m
  logical, allocatable :: do_soag_m(:)

  ! Species presence in modes
  integer, allocatable :: ido_so4a_m(:), ido_nh4a_m(:), ido_soaa_m(:, :)
  integer :: ntot_soamode_m

  ! pcage and pcarbon
  integer :: modetoo_pcage
  integer :: modeptr_pcarbon_m

  ! Mass-to-volume conversion factors
  real(kind_phys) :: fac_m2v_nh4, fac_m2v_so4
  real(kind_phys), allocatable :: fac_m2v_soa(:)

  real(kind_phys), allocatable :: fac_m2v_pcarbon(:)

  ! SOA/POA molecular weights from host model
  real(kind_phys), allocatable :: mw_soa_host_m(:), mw_poa_host_m(:)

  ! Host-provided physical constants:
  real(kind_phys) :: rair_m, mwdry_m
  real(kind_phys) :: rgas_m

contains

  ! initialize gas-aerosol exchange module
  ! store species indices and mode metadata
  ! compute aging/MW conversion factors
  !
  ! Author: R. Easter
  subroutine modal_aero_gasaerexch_init( &
    ntot_amode, nsoa, npoa, nspec_max, &
    nspec_amode, &
    modeptr_pcarbon, modeptr_accum, &
    alnsg_amode, &
    specmw_amode, specdens_amode, spechygro, &
    idx_h2so4, idx_nh3, idx_msa, &
    idx_soag, &
    idx_so4_a, idx_nh4_a, &
    idx_soa_a, idx_pom_a, &
    idx_num, idx_mass, pcnst_in, &
    nspecfrm_pcage_in, &
    lspecfrm_pcage_in, lspectoo_pcage_in, &
    mw_soa_host, mw_poa_host, &
    rair, mwdry, r_universal, &
    errmsg, errflg)

    ! arguments
    integer, intent(in) :: ntot_amode
    integer, intent(in) :: nsoa
    integer, intent(in) :: npoa
    integer, intent(in) :: nspec_max
    integer, intent(in) :: nspec_amode(:)
    integer, intent(in) :: modeptr_pcarbon
    integer, intent(in) :: modeptr_accum
    real(kind_phys), intent(in) :: alnsg_amode(:)
    real(kind_phys), intent(in) :: specmw_amode(:, :)
    real(kind_phys), intent(in) :: specdens_amode(:, :)
    real(kind_phys), intent(in) :: spechygro(:, :)
    integer, intent(in) :: idx_h2so4
    integer, intent(in) :: idx_nh3
    integer, intent(in) :: idx_msa
    integer, intent(in) :: idx_soag(:)
    integer, intent(in) :: idx_so4_a(:)
    integer, intent(in) :: idx_nh4_a(:)
    integer, intent(in) :: idx_soa_a(:, :)
    integer, intent(in) :: idx_pom_a(:, :)
    integer, intent(in) :: idx_num(:)
    integer, intent(in) :: idx_mass(:, :)
    integer, intent(in) :: pcnst_in           ! total number of constituents (for range checks)
    integer, intent(in) :: nspecfrm_pcage_in
    integer, intent(in) :: lspecfrm_pcage_in(:)   ! pcnst-space
    integer, intent(in) :: lspectoo_pcage_in(:)   ! pcnst-space
    real(kind_phys), intent(in) :: mw_soa_host(:)
    real(kind_phys), intent(in) :: mw_poa_host(:)
    real(kind_phys), intent(in) :: rair               ! dry-air gas constant from host (J/K/kg)
    real(kind_phys), intent(in) :: mwdry              ! dry-air molecular weight from host (kg/kmol)
    real(kind_phys), intent(in) :: r_universal        ! universal gas constant from host (J/K/kmol)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! local
    integer :: jsoa, l, l1, l2, n
    real(kind_phys) :: tmp2

    errmsg = ''
    errflg = 0

    ! Store configuration
    ntot_amode_m = ntot_amode
    nsoa_m = nsoa
    npoa_m = npoa
    nspec_max_m = nspec_max

    ! Allocate and store mode configuration arrays
    allocate (nspec_amode_m(ntot_amode))
    nspec_amode_m(:) = nspec_amode(1:ntot_amode)

    allocate (alnsg_amode_m(ntot_amode))
    alnsg_amode_m(:) = alnsg_amode(1:ntot_amode)

    allocate (specmw_amode_m(nspec_max, ntot_amode))
    specmw_amode_m(:, :) = specmw_amode(1:nspec_max, 1:ntot_amode)

    allocate (specdens_amode_m(nspec_max, ntot_amode))
    specdens_amode_m(:, :) = specdens_amode(1:nspec_max, 1:ntot_amode)

    ! Store species indices
    idx_h2so4_m = idx_h2so4
    idx_nh3_m = idx_nh3
    idx_msa_m = idx_msa

    allocate (idx_soag_m(nsoa))
    idx_soag_m(:) = idx_soag(1:nsoa)

    allocate (idx_so4_a_m(ntot_amode))
    idx_so4_a_m(:) = idx_so4_a(1:ntot_amode)

    allocate (idx_nh4_a_m(ntot_amode))
    idx_nh4_a_m(:) = idx_nh4_a(1:ntot_amode)

    allocate (idx_soa_a_m(ntot_amode, nsoa))
    idx_soa_a_m(:, :) = idx_soa_a(1:ntot_amode, 1:nsoa)

    allocate (idx_pom_a_m(ntot_amode, npoa))
    idx_pom_a_m(:, :) = idx_pom_a(1:ntot_amode, 1:npoa)

    allocate (idx_num_m(ntot_amode))
    idx_num_m(:) = idx_num(1:ntot_amode)

    allocate (idx_mass_m(nspec_max, ntot_amode))
    idx_mass_m(:, :) = idx_mass(1:nspec_max, 1:ntot_amode)

    ! Store pcarbon mode pointer
    modeptr_pcarbon_m = modeptr_pcarbon

    ! Store molecular weights
    allocate (mw_soa_host_m(nsoa))
    mw_soa_host_m(:) = mw_soa_host(1:nsoa)

    allocate (mw_poa_host_m(npoa))
    mw_poa_host_m(:) = mw_poa_host(1:npoa)

    ! Store host physical constants
    rair_m = rair
    mwdry_m = mwdry
    rgas_m = r_universal*1.0e-3_kind_phys  ! J/K/kmol -> J/K/mol

    ! Validate H2SO4 index (required species)
    if ((idx_h2so4 <= 0) .or. (idx_h2so4 > pcnst_in)) then
      write (errmsg, '(a,i7)') &
        'modal_aero_gasaerexch_init -- cannot find H2SO4 species, idx=', idx_h2so4
      errflg = 1
      return
    end if

    ! Compute species presence flags
    do_nh4g_m = .false.
    if ((idx_nh3 > 0) .and. (idx_nh3 <= pcnst_in)) do_nh4g_m = .true.

    do_msag_m = .false.
    if ((idx_msa > 0) .and. (idx_msa <= pcnst_in)) do_msag_m = .true.

    allocate (do_soag_m(nsoa))
    do_soag_any_m = .false.
    do_soag_m(:) = .false.
    do jsoa = 1, nsoa
      if ((idx_soag(jsoa) > 0) .and. (idx_soag(jsoa) <= pcnst_in)) then
        do_soag_any_m = .true.
        do_soag_m(jsoa) = .true.
      end if
    end do

    ! Compute ido arrays (species presence in modes)
    allocate (ido_so4a_m(ntot_amode))
    allocate (ido_nh4a_m(ntot_amode))
    allocate (ido_soaa_m(ntot_amode, nsoa))
    ido_so4a_m(:) = 0
    ido_nh4a_m(:) = 0
    ido_soaa_m(:, :) = 0

    ntot_soamode_m = 0
    do n = 1, ntot_amode
      l = idx_so4_a(n)
      if ((l > 0) .and. (l <= pcnst_in)) then
        ido_so4a_m(n) = 1
        if (do_nh4g_m) then
          l = idx_nh4_a(n)
          if ((l > 0) .and. (l <= pcnst_in)) then
            ido_nh4a_m(n) = 1
          end if
        end if
      end if

      do jsoa = 1, nsoa
        if (do_soag_m(jsoa)) then
          l = idx_soa_a(n, jsoa)
          if ((l > 0) .and. (l <= pcnst_in)) then
            ido_soaa_m(n, jsoa) = 1
            ntot_soamode_m = n
          end if
        end if
      end do ! jsoa
    end do ! n

    !
    !   define "from mode" and "to mode" for primary carbon aging
    !
    !   skip (turn off) aging if either is absent,
    !      or if accum mode so4 is absent
    !
    maxspec_pcage = nspec_max
    allocate (lspecfrm_pcage(maxspec_pcage))
    allocate (lspectoo_pcage(maxspec_pcage))
    allocate (soa_equivso4_factor(nsoa))
    allocate (fac_m2v_soa(nsoa))
    allocate (fac_m2v_pcarbon(nspec_max))

    lspecfrm_pcage(:) = 0
    lspectoo_pcage(:) = 0

    modefrm_pcage = -999888777
    modetoo_pcage = -999888777
    nspecfrm_pcage = 0

    if ((modeptr_pcarbon > 0) .and. (modeptr_accum > 0)) then
      l = idx_so4_a(modeptr_accum)
      if ((l > 0) .and. (l <= pcnst_in)) then
        modefrm_pcage = modeptr_pcarbon
        modetoo_pcage = modeptr_accum

        nspecfrm_pcage = nspecfrm_pcage_in
        lspecfrm_pcage(1:nspecfrm_pcage) = lspecfrm_pcage_in(1:nspecfrm_pcage)
        lspectoo_pcage(1:nspecfrm_pcage) = lspectoo_pcage_in(1:nspecfrm_pcage)
      end if
    end if

    if (do_soag_any_m) ntot_soamode_m = max(ntot_soamode_m, modefrm_pcage)

    ! Modify ido arrays for pcage mode
    if (modefrm_pcage > 0) then
      ido_so4a_m(modefrm_pcage) = 2
      if (ido_nh4a_m(modetoo_pcage) == 1) ido_nh4a_m(modefrm_pcage) = 2
      do jsoa = 1, nsoa
        if (ido_soaa_m(modetoo_pcage, jsoa) == 1) ido_soaa_m(modefrm_pcage, jsoa) = 2
      end do
    end if

    ! set for used in aging calcs:
    !    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
    !    soa_equivso4_factor(:)
    soa_equivso4_factor = 0.0_kind_phys
    if (modefrm_pcage > 0) then
      n = modeptr_accum
      l2 = -1
      do l1 = 1, nspec_amode(n)
        if (idx_mass(l1, n) == idx_so4_a(n)) then
!               l2 = lspectype_amode(l1,n)
          l2 = l1
!               fac_m2v_so4 = specmw_amode(l2) / specdens_amode(l2)
          fac_m2v_so4 = specmw_amode(l1, n)/specdens_amode(l1, n)
!               tmp2 = spechygro(l2)
          tmp2 = spechygro(l1, n)

        end if
      end do
      if (l2 <= 0) then
        errmsg = 'modal_aero_gasaerexch_init error a002 finding accum. so4'
        errflg = 1
        return
      end if

      l2 = -1
      if (idx_nh4_a(n) > 0) then
        do l1 = 1, nspec_amode(n)
          if (idx_mass(l1, n) == idx_nh4_a(n)) then
!                  l2 = lspectype_amode(l1,n)
            l2 = l1
!                  fac_m2v_nh4 = specmw_amode(l2) / specdens_amode(l2)
            fac_m2v_nh4 = specmw_amode(l1, n)/specdens_amode(l1, n)

          end if
        end do
        if (l2 <= 0) then
          errmsg = 'modal_aero_gasaerexch_init error a002 finding accum. nh4'
          errflg = 1
          return
        end if
      else
        fac_m2v_nh4 = fac_m2v_so4
      end if

      do jsoa = 1, nsoa
        l2 = -1
        if (idx_soa_a(n, jsoa) <= 0) then
          write (errmsg, '(a,i4)') 'modal_aero_gasaerexch_init error a001 finding accum. jsoa =', jsoa
          errflg = 1
          return
        end if
        do l1 = 1, nspec_amode(n)
          if (idx_mass(l1, n) == idx_soa_a(n, jsoa)) then
!                  l2 = lspectype_amode(l1,n)
            l2 = l1
!                  fac_m2v_soa(jsoa) = specmw_amode(l2) / specdens_amode(l2)
            fac_m2v_soa(jsoa) = specmw_amode(l1, n)/specdens_amode(l1, n)
!                  soa_equivso4_factor(jsoa) = spechygro(l2)/tmp2
            soa_equivso4_factor(jsoa) = spechygro(l1, n)/tmp2
          end if
        end do
        if (l2 <= 0) then
          write (errmsg, '(a,i4)') 'modal_aero_gasaerexch_init error a002 finding accum. jsoa =', jsoa
          errflg = 1
          return
        end if
      end do

      fac_m2v_pcarbon(:) = 0.0_kind_phys
      n = modeptr_pcarbon
      do l = 1, nspec_amode(n)
!            l2 = lspectype_amode(l,n)
!      fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!           [m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!            fac_m2v_pcarbon(l) = specmw_amode(l2) / specdens_amode(l2)
        fac_m2v_pcarbon(l) = specmw_amode(l, n)/specdens_amode(l, n)
      end do
    end if

  end subroutine modal_aero_gasaerexch_init

  subroutine modal_aero_gasaerexch_run( &
    ncol, pver, deltat, top_lev, &
    loffset, &
    aero_props, &
    t, pmid, pdel, gravit, &
    troplev, &
    dgncur_a, dgncur_awet, &
    use_sulfeq, sulfeq, &
    num_q, &
    q, &
    dqdt, dotend, qsrflx_gaexch, &
    dqdt_soa_mode, &
    errmsg, errflg)
    use aerosol_properties_mod, only: aerosol_properties

    integer, intent(in)    :: ncol
    integer, intent(in)    :: pver
    real(kind_phys), intent(in)    :: deltat       ! time step [s]
    integer, intent(in)    :: top_lev              ! top level for aerosol processes
    integer, intent(in)    :: loffset              ! offset to convert pcnst-space to vmr-space [index]
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object (gas uptake rates)
    integer, intent(in)    :: troplev(:)           ! (ncol) tropopause vertical index [index]
    real(kind_phys), intent(in)    :: t(:, :)               ! (ncol,pver) temperature [K]
    real(kind_phys), intent(in)    :: pmid(:, :)            ! (ncol,pver) pressure [Pa]
    real(kind_phys), intent(in)    :: pdel(:, :)            ! (ncol,pver) pressure thickness [Pa]
    real(kind_phys), intent(in)    :: gravit               ! gravitational acceleration [m s-2]
    real(kind_phys), intent(in)    :: dgncur_a(:, :, :)      ! (ncol,pver,ntot_amode) dry diameter
    real(kind_phys), intent(in)    :: dgncur_awet(:, :, :)   ! (ncol,pver,ntot_amode) wet diameter
    logical, intent(in)    :: use_sulfeq           ! whether to use strat equilibrium
    real(kind_phys), intent(in)    :: sulfeq(:, :, :)        ! (ncol,pver,ntot_amode) sulfeq values
    integer, intent(in)    :: num_q                ! # of species in vmr array (= gas_pcnst)
    real(kind_phys), intent(in)    :: q(:, :, :)             ! (ncol,pver,num_q) tracer VMR
    ! *** MUST BE  #/kmol-air for number
    ! *** MUST BE mol/mol-air for mass
    real(kind_phys), intent(out)   :: dqdt(:, :, :)          ! (ncol,pver,num_q) tendencies
    logical, intent(out)   :: dotend(:)            ! (num_q) which species have tendencies
    real(kind_phys), intent(out)   :: qsrflx_gaexch(:, :)   ! (ncol,num_q) column-integrated gas-aerosol
    ! exchange source/sink (kg/m2/s, pre adv_mass/mwdry
    ! scaling) for the _sfgaex1 diagnostic.
    ! Accumulated per-term here because the per-mode
    ! and primary-carbon-aging contributions must be
    ! summed separately to stay bfb with CAM.
    real(kind_phys), intent(out)   :: dqdt_soa_mode(:, :, :, :) ! (ncol,pver,ntot_amode,nsoa) soa aerosol
    ! tendencies of each mode from soa gas uptake
    ! (mol/mol/s), before primary-carbon aging,
    ! for the soa condensation/evaporation diagnostics
    character(len=*), intent(out) :: errmsg
    integer, intent(out)   :: errflg

    ! computes TMR (tracer mixing ratio) tendencies for gas condensation
    !    onto aerosol particles
    !
    ! this version does condensation of H2SO4, NH3, and MSA, both treated as
    ! completely non-volatile (gas --> aerosol, but no aerosol --> gas)
    !    gas H2SO4 goes to aerosol SO4
    !    gas MSA (if present) goes to aerosol SO4
    !       aerosol MSA is not distinguished from aerosol SO4
    !    gas NH3 (if present) goes to aerosol NH4
    !       if gas NH3 is not present, then ????

    ! local variables
    integer, parameter :: method_soa = method_soa_reversible
    !     method_soa=0 is no uptake
    !     method_soa=1 is irreversible uptake done like h2so4 uptake
    !     method_soa=2 is reversible uptake using subr aero_soaexch
    !     (see aero_gasaerexch_core)

    integer :: i, iq
    integer :: ido_so4a(ntot_amode_m), ido_nh4a(ntot_amode_m)
    integer ::  ido_soaa(ntot_amode_m, nsoa_m)
    integer :: j, jsoa
    integer :: k
    integer :: l, lsfrm, lstoo
    integer :: l_so4g, l_nh4g, l_msag
    integer :: l_soag(nsoa_m)
    integer :: n, ntot_soamode

    ! Local offset-adjusted index arrays (pcnst-space - loffset = vmr space)
    integer :: idx_so4_a_q(ntot_amode_m), idx_nh4_a_q(ntot_amode_m)
    integer :: idx_soa_a_q(ntot_amode_m, nsoa_m), idx_pom_a_q(ntot_amode_m, npoa_m)
    integer :: idx_mass_q(nspec_max_m, ntot_amode_m)
    integer :: lspecfrm_q(maxspec_pcage), lspectoo_q(maxspec_pcage)

    real(kind_phys) :: avg_uprt_nh4, avg_uprt_so4
    real(kind_phys) :: deltatxx
    real(kind_phys) :: dqdt_nh4(ntot_amode_m), dqdt_so4(ntot_amode_m)
    real(kind_phys) :: dqdt_soa(ntot_amode_m, nsoa_m)
    real(kind_phys) :: fac_volsfc_pcarbon
    real(kind_phys) :: fgain_nh4(ntot_amode_m), fgain_so4(ntot_amode_m)
    real(kind_phys)  :: mw_poa_host(npoa_m)          ! molec wght of poa used in host code
    real(kind_phys)  :: mw_soa_host(nsoa_m)          ! molec wght of soa used in host code
    real(kind_phys) :: qold_nh4(ntot_amode_m), qold_so4(ntot_amode_m)
    real(kind_phys) :: qold_poa(ntot_amode_m, npoa_m)
    real(kind_phys) :: qold_soa(ntot_amode_m, nsoa_m)
    real(kind_phys) :: qold_soag(nsoa_m)
    real(kind_phys) :: sum_dqdt_msa, sum_dqdt_so4
    real(kind_phys) :: sum_dqdt_soa(nsoa_m)
    real(kind_phys) :: sum_dqdt_nh4, sum_dqdt_nh4_b
    real(kind_phys) :: sum_uprt_nh4, sum_uprt_so4
    real(kind_phys) :: pdel_fac        ! pdel/gravit (kg/m2 per Pa), for column-integrated diagnostics
    real(kind_phys) :: tmp1, tmp2, tmpa
    real(kind_phys) :: uptkrate(ncol, pver, ntot_amode_m)
    real(kind_phys) :: uptkratebb(ntot_amode_m)
    ! gas-to-aerosol mass transfer rates (1/s)
    real(kind_phys) :: vol_core, vol_shell
    real(kind_phys) :: xferfrac_pcage, xferfrac_max
    real(kind_phys) :: xferrate

    logical  :: do_msag         ! true if msa gas is a species
    logical  :: do_nh4g         ! true if nh3 gas is a species
    logical  :: do_soag_any         ! true if soa gas is a species
    logical  :: do_soag(nsoa_m)       ! true if soa gas is a species
    logical  :: do_so4a(ntot_amode_m)       ! true if mode takes up h2so4 gas
    logical  :: do_nh4a(ntot_amode_m)       ! true if mode takes up nh3 gas
    logical  :: do_soamode(ntot_amode_m)    ! true if mode takes part in soa exchange
    logical  :: do_soa(ntot_amode_m, nsoa_m) ! true if soa species takes up soa gas in mode
    logical  :: do_poa(ntot_amode_m)        ! true if poa is included in the soa exchange of a mode

!----------------------------------------------------------------------

    errmsg = ''
    errflg = 0

! set gas species indices from module-level storage, applying -loffset
! to convert pcnst-space to vmr (gas_pcnst) space
    l_so4g = idx_h2so4_m - loffset
    l_nh4g = idx_nh3_m - loffset
    l_msag = idx_msa_m - loffset
    do_nh4g = do_nh4g_m
    do_msag = do_msag_m
    do_soag_any = do_soag_any_m
    do_soag(:) = do_soag_m(:)
    do jsoa = 1, nsoa_m
      l_soag(jsoa) = idx_soag_m(jsoa) - loffset
    end do

! compute offset-adjusted per-mode index arrays
    idx_so4_a_q(:) = idx_so4_a_m(:) - loffset
    idx_nh4_a_q(:) = idx_nh4_a_m(:) - loffset
    idx_soa_a_q(:, :) = idx_soa_a_m(:, :) - loffset
    idx_pom_a_q(:, :) = idx_pom_a_m(:, :) - loffset
    idx_mass_q(:, :) = idx_mass_m(:, :) - loffset
    do iq = 1, nspecfrm_pcage
      lspecfrm_q(iq) = lspecfrm_pcage(iq) - loffset
      lspectoo_q(iq) = lspectoo_pcage(iq)
      if (lspectoo_q(iq) > 0) lspectoo_q(iq) = lspectoo_q(iq) - loffset
    end do

! copy ido arrays from module-level storage
    ido_so4a(:) = ido_so4a_m(:)
    ido_nh4a(:) = ido_nh4a_m(:)
    ido_soaa(:, :) = ido_soaa_m(:, :)
    ntot_soamode = ntot_soamode_m

! set molecular weights from module-level storage
    mw_soa_host(:) = mw_soa_host_m(:)
    mw_poa_host(:) = mw_poa_host_m(:)

! set tendency flags
    dotend(:) = .false.

    dotend(l_so4g) = .true.
    if (do_nh4g) dotend(l_nh4g) = .true.
    if (do_msag) dotend(l_msag) = .true.
    do jsoa = 1, nsoa_m
      if (do_soag(jsoa)) dotend(l_soag(jsoa)) = .true.
    end do

    do n = 1, ntot_amode_m
      if (ido_so4a(n) == 1) then
        l = idx_so4_a_q(n)
        dotend(l) = .true.
        if (do_nh4g) then
          if (ido_nh4a(n) == 1) then
            l = idx_nh4_a_q(n)
            dotend(l) = .true.
          end if
        end if
      end if

      do jsoa = 1, nsoa_m
        if (do_soag(jsoa)) then
          if (ido_soaa(n, jsoa) == 1) then
            l = idx_soa_a_q(n, jsoa)
            dotend(l) = .true.
          end if
        end if
      end do ! jsoa
    end do ! n

    if (modefrm_pcage > 0) then
      do iq = 1, nspecfrm_pcage
        lsfrm = lspecfrm_q(iq)
        lstoo = lspectoo_q(iq)
        if ((lsfrm > 0) .and. (lsfrm <= num_q)) then
          dotend(lsfrm) = .true.
          if ((lstoo > 0) .and. (lstoo <= num_q)) then
            dotend(lstoo) = .true.
          end if
        end if
      end do

      n = modeptr_pcarbon_m
      fac_volsfc_pcarbon = exp(2.5_kind_phys*(alnsg_amode_m(n)**2))
      xferfrac_max = 1.0_kind_phys - 10.0_kind_phys*epsilon(1.0_kind_phys)   ! 1-eps
    end if

! modes taking up h2so4 and nh3 gas
    do n = 1, ntot_amode_m
      do_so4a(n) = (ido_so4a(n) > 0)
      do_nh4a(n) = (ido_nh4a(n) > 0)
    end do

! modes and species taking up soa gas
    do n = 1, ntot_amode_m
      do_soamode(n) = (n <= ntot_soamode)
      do jsoa = 1, nsoa_m
        do_soa(n, jsoa) = (ido_soaa(n, jsoa) > 0)
      end do
      ! in previous code versions with nsoa=1,
      !    qold_poa was non-zero (i.e., loaded from q) only when ido_soaa(n)=1
      ! thus qold_poa=0 for the primary carbon mode which has ido_soaa=2
      ! this is probably not how it should be
      do_poa(n) = any(ido_soaa(n, :) == 1)
    end do

! zero out tendencies
    dqdt(:, :, :) = 0.0_kind_phys
    qsrflx_gaexch(:, :) = 0.0_kind_phys
    dqdt_soa_mode(:, :, :, :) = 0.0_kind_phys

! compute gas-to-aerosol mass transfer rates
    call gas_aer_uptkrates(ncol, pver, top_lev, &
                           loffset, aero_props, &
                           q, t, pmid, &
                           dgncur_awet, uptkrate)

! use this for tendency calcs to avoid generating very small negative values
    deltatxx = deltat*(1.0_kind_phys + 1.0e-15_kind_phys)

    do k = top_lev, pver
      do i = 1, ncol

        do n = 1, ntot_amode_m
          uptkratebb(n) = uptkrate(i, k, n)
          if (ido_so4a(n) == 1) then
            qold_so4(n) = q(i, k, idx_so4_a_q(n))
          else
            qold_so4(n) = 0.0_kind_phys
          end if
          if (ido_nh4a(n) == 1) then
            qold_nh4(n) = q(i, k, idx_nh4_a_q(n))
          else
            qold_nh4(n) = 0.0_kind_phys
          end if
        end do ! n

!   fgain_so4(n) = fraction of total h2so4 uptake going to mode n
!   fgain_nh4(n) = fraction of total  nh3  uptake going to mode n
!   avg_uprt_so4/nh4 = uptake amount (fraction of gas uptaken) over deltat
!   (nh3_diffus_ratio is for gas diffusivity (nh3/h2so4);
!    differences in fuch-sutugin and accom coef ignored)
        call aero_gasaerexch_partition(ntot_amode_m, do_so4a, uptkratebb, h2so4_diffus_ratio, deltatxx, &
                                       fgain_so4, sum_uprt_so4, avg_uprt_so4)
        call aero_gasaerexch_partition(ntot_amode_m, do_nh4a, uptkratebb, nh3_diffus_ratio, deltatxx, &
                                       fgain_nh4, sum_uprt_nh4, avg_uprt_nh4)

!   sum_dqdt_so4 = so4_a tendency from h2so4 gas uptake (mol/mol/s)
!   sum_dqdt_msa = msa_a tendency from msa   gas uptake (mol/mol/s)
!   sum_dqdt_nh4 = nh4_a tendency from nh3   gas uptake (mol/mol/s)
        sum_dqdt_so4 = q(i, k, l_so4g)*avg_uprt_so4
        if (do_msag) then
          sum_dqdt_msa = q(i, k, l_msag)*avg_uprt_so4
        else
          sum_dqdt_msa = 0.0_kind_phys
        end if
        if (do_nh4g) then
          sum_dqdt_nh4 = q(i, k, l_nh4g)*avg_uprt_nh4
        else
          sum_dqdt_nh4 = 0.0_kind_phys
        end if

        if (use_sulfeq .and. (k <= troplev(i))) then
          !   compute TMR tendencies for so4 interstial aerosol due to reversible gas uptake
          !   only above the tropopause
          !   (dqdt_so4 is left unchanged in modes that do not take up h2so4)
          call aero_gasaerexch_reversible(ntot_amode_m, do_so4a, uptkratebb, sulfeq(i, k, :), deltatxx, &
                                          sum_uprt_so4, q(i, k, l_so4g), qold_so4, &
                                          dqdt_so4, sum_dqdt_so4)
          ! do not allow msa condensation in stratosphere
          ! ( Note that the code for msa has never been used.
          !   The plan was to simulate msa(g), treat it as non-volatile (like h2so4(g)),
          !   and treat condensed msa as sulfate, so just one additional tracer. )
          if (do_msag) sum_dqdt_msa = 0.0_kind_phys

        else
          !   compute TMR tendencies for so4 interstial aerosol due to simple gas uptake
          do n = 1, ntot_amode_m
            dqdt_so4(n) = fgain_so4(n)*(sum_dqdt_so4 + sum_dqdt_msa)
          end do
        end if

        !   compute TMR tendencies for nh4 interstial aerosol due to simple gas uptake
        !   but force nh4/so4 molar ratio <= 2
        sum_dqdt_nh4_b = 0.0_kind_phys
        dqdt_nh4(:) = 0._kind_phys
        if (do_nh4g) then
          call aero_gasaerexch_nh4_limit(ntot_amode_m, fgain_nh4, sum_dqdt_nh4, deltat, deltatxx, &
                                         qold_nh4, qold_so4, dqdt_so4, dqdt_nh4, sum_dqdt_nh4_b)
        end if

        if (do_soag_any) then
          !   compute TMR tendencies for soag and soa interstial aerosol
          !   (soa loaded only for species with ido_soaa = 1; poa only in modes with do_poa)
          do jsoa = 1, nsoa_m
            if (do_soag(jsoa)) then
              qold_soag(jsoa) = q(i, k, l_soag(jsoa))
            else
              qold_soag(jsoa) = 0.0_kind_phys
            end if
          end do
          do n = 1, ntot_amode_m
            do jsoa = 1, nsoa_m
              if (ido_soaa(n, jsoa) == 1) then
                qold_soa(n, jsoa) = q(i, k, idx_soa_a_q(n, jsoa))
              else
                qold_soa(n, jsoa) = 0.0_kind_phys
              end if
            end do
            do j = 1, npoa_m
              l = idx_pom_a_q(n, j)
              if (do_poa(n) .and. (l > 0)) then
                qold_poa(n, j) = q(i, k, l)
              else
                qold_poa(n, j) = 0.0_kind_phys
              end if
            end do
          end do

          call aero_gasaerexch_soa(method_soa, deltat, t(i, k), pmid(i, k), rgas_m, &
                                   ntot_amode_m, nsoa_m, npoa_m, do_soamode, do_soa, do_soag, uptkratebb, &
                                   mw_poa_host, mw_soa_host, &
                                   qold_soag, qold_soa, qold_poa, &
                                   dqdt_soa, sum_dqdt_soa, errmsg, errflg)
          if (errflg /= 0) return
        else
          dqdt_soa(:, :) = 0.0_kind_phys
          sum_dqdt_soa(:) = 0.0_kind_phys
        end if
        dqdt_soa_mode(i, k, :, :) = dqdt_soa(:, :)

        pdel_fac = pdel(i, k)/gravit
        do n = 1, ntot_amode_m
          if (ido_so4a(n) == 1) then
            l = idx_so4_a_q(n)
            dqdt(i, k, l) = dqdt_so4(n)
            qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt_so4(n)*pdel_fac
          end if

          if (do_nh4g) then
            if (ido_nh4a(n) == 1) then
              l = idx_nh4_a_q(n)
              dqdt(i, k, l) = dqdt_nh4(n)
              qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt_nh4(n)*pdel_fac
            end if
          end if

          do jsoa = 1, nsoa_m
            if (do_soag(jsoa)) then
              if (ido_soaa(n, jsoa) == 1) then
                l = idx_soa_a_q(n, jsoa)
                dqdt(i, k, l) = dqdt_soa(n, jsoa) !calculated by  aero_soaexch for method_soa=2
                qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt_soa(n, jsoa)*pdel_fac
              end if
            end if
          end do
        end do ! n

        !   compute TMR tendencies for h2so4, nh3, and msa gas
        !   due to simple gas uptake
        l = l_so4g
        dqdt(i, k, l) = -sum_dqdt_so4
        qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt(i, k, l)*pdel_fac

        if (do_msag) then
          l = l_msag
          dqdt(i, k, l) = -sum_dqdt_msa
          qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt(i, k, l)*pdel_fac
        end if

        if (do_nh4g) then
          l = l_nh4g
          dqdt(i, k, l) = -sum_dqdt_nh4_b
          qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt(i, k, l)*pdel_fac
        end if

        do jsoa = 1, nsoa_m
          if (do_soag(jsoa)) then
            l = l_soag(jsoa)
            dqdt(i, k, l) = -sum_dqdt_soa(jsoa)
! dqdt for gas is negative of the sum of dqdt for aerosol soa species in each mode: Manish
            qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt(i, k, l)*pdel_fac
          end if
        end do

!   compute TMR tendencies associated with primary carbon aging
        if (modefrm_pcage > 0) then
          n = modeptr_pcarbon_m
          tmpa = 0.0_kind_phys
          do jsoa = 1, nsoa_m
            tmpa = tmpa + dqdt_soa(n, jsoa)*fac_m2v_soa(jsoa)*soa_equivso4_factor(jsoa)
          end do
          vol_shell = deltat* &
                      (dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 + tmpa)
          vol_core = 0.0_kind_phys
          do l = 1, nspec_amode_m(n)
            vol_core = vol_core + &
                       q(i, k, idx_mass_q(l, n))*fac_m2v_pcarbon(l)
          end do

          !   ratio1 = vol_shell/vol_core =
          !      actual hygroscopic-shell-volume/carbon-core-volume after gas uptake
          !   ratio2 = 6.0_kind_phys*dr_so4_monolayers_pcage/(dgncur_a*fac_volsfc_pcarbon)
          !      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume
          !      The 6.0/(dgncur_a*fac_volsfc_pcarbon) = (mode-surface-area/mode-volume)
          !   Note that vol_shell includes both so4+nh4 AND soa as "equivalent so4",
          !      The soa_equivso4_factor accounts for the lower hygroscopicity of soa.
          !
          !   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
          !   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow

          tmp1 = vol_shell*dgncur_a(i, k, n)*fac_volsfc_pcarbon
          tmp2 = max(6.0_kind_phys*dr_so4_monolayers_pcage*vol_core, 0.0_kind_phys)
          if (tmp1 >= tmp2) then
            xferfrac_pcage = xferfrac_max
          else
            xferfrac_pcage = min(tmp1/tmp2, xferfrac_max)
          end if

          if (xferfrac_pcage > 0.0_kind_phys) then
            do iq = 1, nspecfrm_pcage
              lsfrm = lspecfrm_q(iq)
              lstoo = lspectoo_q(iq)
              xferrate = (xferfrac_pcage/deltat)*q(i, k, lsfrm)
              dqdt(i, k, lsfrm) = dqdt(i, k, lsfrm) - xferrate
              qsrflx_gaexch(i, lsfrm) = qsrflx_gaexch(i, lsfrm) - xferrate*pdel_fac
              if ((lstoo > 0) .and. (lstoo <= num_q)) then
                dqdt(i, k, lstoo) = dqdt(i, k, lstoo) + xferrate
                qsrflx_gaexch(i, lstoo) = qsrflx_gaexch(i, lstoo) + xferrate*pdel_fac
              end if
            end do

            if (ido_so4a(modetoo_pcage) > 0) then
              l = idx_so4_a_q(modetoo_pcage)
              dqdt(i, k, l) = dqdt(i, k, l) + dqdt_so4(modefrm_pcage)
              qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt_so4(modefrm_pcage)*pdel_fac
            end if

            if (ido_nh4a(modetoo_pcage) > 0) then
              l = idx_nh4_a_q(modetoo_pcage)
              dqdt(i, k, l) = dqdt(i, k, l) + dqdt_nh4(modefrm_pcage)
              qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt_nh4(modefrm_pcage)*pdel_fac
            end if

            do jsoa = 1, nsoa_m
              if (ido_soaa(modetoo_pcage, jsoa) > 0) then
                l = idx_soa_a_q(modetoo_pcage, jsoa)
                dqdt(i, k, l) = dqdt(i, k, l) + dqdt_soa(modefrm_pcage, jsoa)
                qsrflx_gaexch(i, l) = qsrflx_gaexch(i, l) + dqdt_soa(modefrm_pcage, jsoa)*pdel_fac
              end if
            end do

          end if

        end if

      end do   ! "i = 1, ncol"
    end do     ! "k = top_lev, pver"

  end subroutine modal_aero_gasaerexch_run

  ! computes gas-to-aerosol mass transfer rates (1/s) of all modes
  ! (see aerosol_properties%gas_uptkrate)
  subroutine gas_aer_uptkrates(ncol, pver, top_lev, &
                               loffset, aero_props, &
                               q, t, pmid, &
                               dgncur_awet, uptkrate)
    use aerosol_properties_mod, only: aerosol_properties

    integer, intent(in) :: ncol                 ! number of atmospheric column
    integer, intent(in) :: pver                 ! number of vertical levels
    integer, intent(in) :: top_lev              ! top level for aerosol processes
    integer, intent(in) :: loffset              ! offset to convert pcnst-space to vmr space
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(kind_phys), intent(in) :: q(:, :, :)             ! (ncol,pver,num_q) Tracer array (mol,#/mol-air)
    real(kind_phys), intent(in) :: t(:, :)               ! (ncol,pver) Temperature in Kelvin
    real(kind_phys), intent(in) :: pmid(:, :)            ! (ncol,pver) Air pressure in Pa
    real(kind_phys), intent(in) :: dgncur_awet(:, :, :)   ! (ncol,pver,ntot_amode_m)

    real(kind_phys), intent(out) :: uptkrate(:, :, :)     ! (ncol,pver,ntot_amode_m)
    ! gas-to-aerosol mass transfer rates (1/s)

! local
    integer :: i, k, n

    real(kind_phys) :: aircon
    real(kind_phys) :: numdens(ncol, pver)       ! mode number density (#/m3)
    real(kind_phys) :: rhoair

    ! 22-aug-2007 rc easter - get number from q array rather
    !    than computing a "bounded" number conc.

    ! outermost loop over all modes
    do n = 1, ntot_amode_m

      numdens(:, :) = 0.0_kind_phys
      do k = top_lev, pver
        do i = 1, ncol
          rhoair = pmid(i, k)/(rair_m*t(i, k))   ! (kg-air/m3)
          !   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air)
          !   so need aircon in (kmol-air/m3)
          aircon = rhoair/mwdry_m            ! (kmol-air/m3)
          numdens(i, k) = q(i, k, idx_num_m(n) - loffset)*aircon
        end do
      end do

      call aero_props%gas_uptkrate(n, ncol, pver, top_lev, numdens, dgncur_awet(:, :, n), &
                                   t, pmid, uptkrate(:, :, n))

    end do

  end subroutine gas_aer_uptkrates

end module modal_aero_gasaerexch
