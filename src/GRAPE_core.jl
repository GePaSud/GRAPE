# GRAPE = General Relativity Accelerometer-based Propagation Environment
# coding in Julia by Jean-Pierre Barriot, July 4, 2025,
# DO NOT DISTRIBUTE WITHOUT PERMISSION TO DO SO
# please do NOT distribute without consent of authors, this is NOT yet open source
# submission to software-X journal in progress
# Windows 11 version, Julia 1.11.5
#------------------------------------------------------------------------------------------------------------------------------------------------
# integration of the motion of a spacecraft (S/C) in a general relativity framework with associated attitude by Fermi-Walker transport,
# in presence (or not) of non-gravitational forces
# non-gravitational forces are either coded wrt to chart coordinates, as the gradient of the energy-momentum tensor, or as locally modeled or measured 
# acceklerations in the S/c bus frame (so the name Accelerometer-based)
#------------------------------------------------------------------------------------------------------------------------------------------------
# symplectic integrator with conservation of the norm of the 4-velocity
# Gauss 5-stages symplectic integrator (Butcher)
# symplecticity is enforced by integrating the equations of motion of the S/C wrt to an instantaneously-at-rest tetrad along the worldline of the S/C
# named Cholesky tetrad in the Julia code
# see O'Leary and Barriot, 2021, https://doi.org/10.1007/s10569-021-10051-7, for the mathematics
# and Barriot, O'Leary and Yan, Journal of Physics, conference Series, 2025, in press, for the exact algorithm
#------------------------------------------------------------------------------------------------------------------------------------------------
# three different times (clocks) are considered for the integration of the equations of motion of the S/C:
# 1/ the coordinate time (also called ephemeris-time), i.e. the time of a non-rotating, not-moving clock outside the gravity well of the Solar System,
# 2/ the proper time of the S/C, the time of a clock co-moving with the S/C
# 3/ an integration time, that is either 1/ or 2/ clocks, or a end-user defined clock, that must be linked by the end-user to either 1/ or 2/ clocks
#------------------------------------------------------------------------------------------------------------------------------------------------
# the code is organized around a set of 26 partially coupled first-order linear equations, that are divided in 3 groups:
# 1/ the first eight equations (1-8 in Julia numbering) described the motion (chart coordinates, wrt to integration time) of the center-of-gravity of a S/C,
# the first equation of motion gives the coordinate time wrt to integration time
# the S/C velocities outputted are 4-velocities (i.e. wrt to proper time) in chart coordinates, but their evolution is also given wrt to integration time
# 2/ equation (9) is the proper time wrt integration time, equation (10) directly gives access to the integration time
# 3/ equations (11-26) are relative to the Fermi-Walker transported co-moving tetrad, that permits the orientation of the S/C wrt distant stars
# one of the 4-vectors of this transported tetrad must be the 4-velocity (so the name co-moving),
# in order to keep the orientation of this terad constant wrt distant stars
# equations (11-26) are not necessary for the integration of equations(1-10), an can be omitted for ephemeris applications
#------------------------------------------------------------------------------------------------------------------------------------------------
# in this code, as an example, we consider the orbit of the Parker Solar probe probe in the Kerr's polar metric
# (see routine Kerr_Metric_Polar) in spatial polar coordinates, (+---) Lorentzian signature,
# that can be easily adapted to any asymptotically flat metric, including a metric defined wrt to Solar System Ephemeris like JPL-DE440
# as the connection connections are computed numerically,
# with the restriction that the first coordinate of integration (and the metric) must be the time coordinate.
# The integration time can be either proper time or coordinate time, or a function of them (see routine evpar2tau_external)
#------------------------------------------------------------------------------------------------------------------------------------------------
# all floating points quantities are coded with **extended precision** using the BigFloat format, do **NOT** suppress this feature,
# it is absolutely needed for a proper modeling of the relativistic effects, see keyword
# setprecision(BigFloat, xxx) # xxx = 53 for 64 bits IEEE754, 113 for 128 bits, 237 for 256 bits IEEE754,
# by using the GNU-MPFR library built natively in Julia
# this is *why* this code was developped in Julia, where BigFloats are, again, natively implemented
# porting this code to other computer languages must be done with extreme care about the intrinsic precision of floating variables on these languages
#------------------------------------------------------------------------------------------------------------------------------------------------
# the code also outputs a pseudo-Doppler observable (as a display), that corresponds to the ratio between
# the derivative of the coordinate time versus integration time divided by the derivative of the proper time versus the integration time,
# this ratio is the second component of the integrated state vector (26 components)
# it corresponds to a Doppler counter onboard the S/C receiving an ingoing wave emitted at infinity
# from a non-rotating source at rest, with no preferred direction, and wrt
# the frequency given in the variable "transmit_frequency_Doppler"
#------------------------------------------------------------------------------------------------------------------------------------------------
# a little bit of windows 11 targeted routines, for computing environment purposes, can be adapted to Linux/Mac without too much trouble
# in some computers, implementing PythonPlot library can be troublesome if done directly from Julia, implement miniconda first
# tests were done with Julia version 1.11.5 
#----------------------------------------------------------- code starts here------------------------------------------------------------------------------
# ------------------------------------------------- first some comfort routines ---------------------------------------------------------------------------
function Call_System(command::String)
# to call the underlying Windows shell
# windows OS only, most of the basic commands should be "cmd /c ....", /c is to exit shell
# command = "cmd /c cd $working_directory"
# system_message = Call_System(command)
# println("system_message= ", system_message)
        if !isascii(command)
           error("input of Call_System must be an all-ascii string")
        end
        split_string = split(command)
        return_string = run(`$split_string`)
    return return_string
end
function prt2(x...)
# Python-like print, with spaces printed between values for ease of use
    N = length(x)
    for i = 1:N-1
        print(x[i], " ")
    end
    println(x[N])
    #write(runlogfile, timetag*"\n")
end
function prt3(x...)
# combines prt3 and prtstr and save in the runlogfile
    N = length(x)
    str = ""
    for i = 1:N-1
        print(x[i], " ")
        str = str*string(x[i])*" "
    end
    println(x[N])
    str = str*string(x[N])*"\n"
    write(runlogfile, str)
    return
end
function Give_Function_Name(Function_Name::Function)
    return Function_Name
end
# ------------------------------------------------- start of physics routines --------------------------------------------------------
function Minkowski_Metric_Cartesian(y::Vector{BigFloat})
# in cartesian local coordinates, NOT in polar form
# y must be provided but is not used, to avoid a possible compilation error (empty entry list)
# Minkowski covariant metric tensor (c^2, -1, -1, -1)
        mt_m = fill(BF0, (4, 4))
        mt_m[1, 1] = c_light^2
        mt_m[2, 2] = -BF1 ; mt_m[3, 3] = -BF1 ; mt_m[4, 4] = -BF1
        # Minkowski contravariant metric tensor
        mt_mi = fill(BF0, (4, 4))
        mt_mi[1, 1] = BF1/mt_m[1, 1] ; mt_mi[2, 2] = BF1/mt_m[2, 2]
        mt_mi[3, 3] = BF1/mt_m[3, 3] ; mt_mi[4, 4] = BF1/mt_m[4, 4]
        return (mt_m, mt_mi)
end
function Minkowski_Metric_Polar(y::Vector{BigFloat})
# for tests, polar coordinates
# t = y[1] = coordinate time (clock at infinity) ; r = y[3] = pseudo-radius
# phi = y[5] = pseudo-longitude ; theta = y[7] = pseudo-colatitude,
        r = y[3] ; theta = y[7]
        sn = sin(theta)
        # Minkowski covariant metric tensor
        mt_s = fill(BF0, (4, 4))
        mt_s[1, 1] = +c_light^2 ; mt_s[2, 2] = -BF1          
        mt_s[3, 3] = -(r*sn)^2    ; mt_s[4, 4] = -r^2
        # Minkowski contravariant metric tensor
        mt_si = fill(BF0, (4, 4))
        mt_si[1, 1] = BF1/mt_s[1, 1] ; mt_si[2, 2] = BF1/mt_s[2, 2]
        mt_si[3, 3] = BF1/mt_s[3, 3] ; mt_si[4, 4] = BF1/mt_s[4, 4]
        return (mt_s, mt_si)
end
function Schwarzschild_Metric_Polar(y::Vector{BigFloat})
# for tests, polar coordinates
# Schwarszchild metric, reduces to previous Minkowski_Metric_Cartesian if mass of star is zero (value rs = 0)
# t = y[1] = coordinate time (clock at infinity) ; r = y[3] = pseudo-radius
# phi = y[5] = pseudo-longitude ; theta = y[7] = pseudo-colatitude,
    r = y[3] ; theta = y[7]
    sn = sin(theta)
    w = BF1-rs/r # rs = Schwarzschild radius
    # Schwarszchild covariant metric tensor
    mt_s = fill(BF0, (4, 4))
    mt_s[1, 1] = +c_light^2*w ; mt_s[2, 2] = -BF1/w          
    mt_s[3, 3] = -(r*sn)^2    ; mt_s[4, 4] = -r^2
    # Schwarszchild contravariant metric tensor
    mt_si = fill(BF0, (4, 4))
    mt_si[1, 1] = BF1/mt_s[1, 1] ; mt_si[2, 2] = BF1/mt_s[2, 2]
    mt_si[3, 3] = BF1/mt_s[3, 3] ; mt_si[4, 4] = BF1/mt_s[4, 4]
    return (mt_s, mt_si)
end
function Schwarzschild_BG_Polar(y::Vector{BigFloat})
# Christoffel symbols for Schwarszchild metric in polar coordinates, reduces to previous Minkowski_Metric_Cartesian if mass of star is zero (value rs = 0)
# t = y[1] = coordinate time (clock at infinity) ; r = y[3] = pseudo-radius
# phi = y[5] = pseudo-longitude ; theta = y[7] = pseudo-colatitude,
# rs = Schwarzschild radius
    BG = fill(BF0, (4, 4, 4))
# indexes t ==1 ; r==2 ; 3 == longitude phi ; 4 == colatitude theta
    r = y[3] ; theta = y[7]
    BG[1, 1, 2] = BG[1, 2, 1] = rs/BF2/r/(r-rs)
    BG[2, 1, 1] = c_light^2*rs*(r-rs)/BF2/r^3
    BG[2, 2, 2] = -rs/BF2/r/(r-rs)
    BG[2, 3, 3] = -(r-rs)*(sin(theta))^2
    BG[2, 4, 4] = -(r-rs)
    BG[3, 2, 3] = BG[3, 3, 2] = BF1/r
    BG[3, 3, 4] = BG[3, 4, 3] = cot(theta)
    BG[4, 3, 3] = -sin(theta)*cos(theta)
    BG[4, 2, 4] = BG[4, 4, 2] = BF1/r
    return BG
end
function Schwarzschild_tetrad(y::Vector{BigFloat})
    (mt_s, mt_si) = Schwarzschild_Metric_Polar(y)
    tetrad = fill(BF0, (4, 4))
    tetradinv = fill(BF0, (4, 4))
    r = y[3] ; theta = y[7]
    sn = sin(theta)
    w = BF1-rs/r # rs = Schwarzschild radius
    tetradinv[1, 1] = sqrt(w) ; tetradinv[2, 2] = BF1/sqrt(w )         
    tetradinv[3, 3] = r*sn    ; tetradinv[4, 4] = r
    tetrad = inv(tetradinv)
    return(tetrad, tetradinv)
end
function TDB_FB(TT::BigFloat)
# TDB (coordinate, or ephemeris time) as a function of TT (Terrestrial Time)
# for tests, not used in this run example
# Kopeikin et al. book (DOI:10.1002/9783527634569), page 715,
# abridged expression from Fairhead and Bretagnon, 1990
TT = Float64(TT)
T = (TT-2451545.0)/36525.0 # TT (Terrestrial Time) in JD
TDB = TT+0.001657  *sin( 628.3076*T+6.2401)
         +0.000022  *sin( 575.3385*T+4.2970)
         +0.000014  *sin(1256.6152*T+6.1969)
         +0.000005  *sin( 606.9777*T+4.0212)
         +0.000005  *sin(  52.9691*T+0.4444)
         +0.000002  *sin(  21.3299*T+5.5431)
         +0.000010*T*sin( 628.3076*T+4.2490)
TDB = parse(BigFloat, string(TDB))
    return TDB
end
function verify_harmonic_coord(y::Vector{BigFloat}, Given_Metric::Function)
# verify if he coordinates are harmonic wrt the metric tensor
# for tests, not used in this run example
# BigGamma # first indice = upper position of the Christoffel symbol,
# the two lower indices are symetrical for holonomic coordinates
   (mt_c, mt_ci) = Given_Metric(y)
   BigGamma = Numerical_Christoffel_Symbols(order_parder_metric, y, Given_Metric)
   verif_harm = fill(BF0, 4)
    for i = 1:4
        for j = 1:4
            for k = 1:4
                verif_harm[i] = verif_harm[i]+mt_ci[j, k]*BigGamma[i, j, k]
            end
        end
    end
    return verif_harm
end
#=
# --------------------for possible use with Newtonian_Metric_Polar-----------------------------------
function init_de430()
# initialisations to SPICE functions (DE430 ephemeris of the Solar System that was used to model PSP orbit),
# see: https://naif.jpl.nasa.gov/naif/toolkit.html and https://pds.nasa.gov/
# for all calls to furnsh, you need to provide the exact path where these files are located
# using SPICE must be active in the main code
   furnsh("naif0012.tls") # leap seconds for UTC to ET
   furnsh("de430.bsp") # ephemeris of main bodies of the Solar System
   furnsh("ast343de430.bsp") # ephemeris of minor bodies
   furnsh("spp_nom_20180812_20251001_v041_RO8.bsp") # reconstructed orbit of PSP (wrt to DE430)
# GM of all planet/moons systems in the Solar System (for use with de430.bsp only)
   BODY_GM = fill(BF0, 10) # units: KM**3/SEC**2
   BODY_GM[1] = parse(BigFloat, "22031.780000") # Mercury
   BODY_GM[2] = parse(BigFloat, "324858.592000") # Venus
   BODY_GM[3] = parse(BigFloat, "403503.235502") # Earth / Moon
   BODY_GM[4] = parse(BigFloat, "42828.375214") # Mars + moons
   BODY_GM[5] = parse(BigFloat, "126712764.800000") # Jupiter + moons
   BODY_GM[6] = parse(BigFloat, "37940585.200000") # Saturn + moons
   BODY_GM[7] = parse(BigFloat, "5794548.600000") # Uranus + moons
   BODY_GM[8] = parse(BigFloat, "6836527.100580") # Neptune + moons
   BODY_GM[9] = parse(BigFloat, "977.000000") # Pluto / Charon
   BODY_GM[10] = parse(BigFloat, "132712440041.939400") # Sun
#--------------------------------------------------------------------------------------
return BODY_GM
end
=#
function Newtonian_Metric_Polar(y::Vector{BigFloat})
    # this is the standard approximation of the metric by Newtonian potential, valid up to 1/c_light^4 order, with no frame-dragging
    r = y[3] ; theta = y[7]
    sn = sin(theta)
    newtonian_potential = mu/r # with a "+" sign, can be replaced by any valid potential model of the Solar System , for exmple derived from DE430
    # see commented routine init_de430()
    mt_s = fill(BF0, (4, 4))
    mt_s[1, 1] = (c_light^2-BF2*newtonian_potential) ; mt_s[2, 2] = -(BF1+BF2/c_light^2*newtonian_potential)        
    mt_s[3, 3] = -(r*sn)^2    ; mt_s[4, 4] = -r^2
    # Newtonian contravariant metric tensor
    mt_si = fill(BF0, (4, 4))
    mt_si[1, 1] = BF1/mt_s[1, 1] ; mt_si[2, 2] = BF1/mt_s[2, 2]
    mt_si[3, 3] = BF1/mt_s[3, 3] ; mt_si[4, 4] = BF1/mt_s[4, 4]
    return (mt_s, mt_si)
# more complicated weak-field metrics, including PPN parameters,
# can be constructed here (see for example the book of Brumberg, https://doi.org/10.1201/9780203756591)
# and using the SPICE kernels: see https://naif.jpl.nasa.gov/pub/naif/misc/toolkit_docs_N0067/C/req/naif_ids.html
end
#-------------------------------this is the routine currently implemented in this code-----------------------------------------------
function Kerr_Metric_Polar(y::Vector{BigFloat})
# Kerr metric, reduces to Schwarzschild_Metric_Polar if the star is not rotating (value as = 0)
# Boyer–Lindquist polar coordinates (technically oblate spheroidal coordinates),
# see catalog of spacetimes, https://doi.org/10.48550/arXiv.0904.4184
# the "as" parameter is the scaled angular momentum, in units of length
# as = Jstar/Mstar/c_light
    r = y[3] ; theta = y[7]
    sn = sin(theta) ; cs = cos(theta)
    BigSigma = r^2+as^2*cs^2
    BigDelta = r^2-rs*r+as^2
    mt_k = fill(BF0, (4, 4))
    mt_k[1, 1] = +c_light^2*(BF1-rs*r/BigSigma)
    mt_k[1, 3] = mt_k[3, 1] = rs*as*r*sn^2*c_light/BigSigma # multiplied by 2 by symmetry
    mt_k[2, 2] = -BigSigma/BigDelta
    mt_k[3, 3] = -(r^2+as^2+rs*as^2*r*sn^2/BigSigma)*sn^2
    mt_k[4, 4] = -BigSigma
    mt_ki = inv(mt_k)
    return (mt_k, mt_ki)
    end
#function Any_SpaceTime_Metric(y)
# you can provide here the metric you want, including metrics build on DE430 JPL ephemeris or others, provided that the first coordinate
# MUST be the coordinate time (TBD, ephemeris_time, i.e. a clock at spatial infinity). This is because coordinate time is the 
# the time that links all events in the Solar System on a common clock
# signature must be (+, -, -, -)
#end
#
function Minkowski_norm2(vector::Vector{BigFloat})
# returns the Minkowskian (Cartesian) squared norm of a 4-vector
    norm2 = c_light^2*vector[1]^2-vector[2]^2-vector[3]^2-vector[4]^2
    return norm2
end
function evpar2tau_external(time_option::Integer, evpar::BigFloat, y::Vector{BigFloat}, Given_Metric::Function)
# gives the proper time as a function imposed by an a priori relationship, EXTERNAL to the code
# for example, the proper time of the spacecraft wrt to an Earth clock
# here tau = sigma, sigma being noted evpar (EVolution_PARameter) as variable in the julia code
# nor y nor Given_Metric are used in this routine as on April 28, 2025, reserved for future use
# a practical example linking terrestrial time (TT) with coordinate time TDB (or ET) can be found in Kopeikin et al., , eq. 9.6, page 715
# at DOI:10.1002/9783527634569
    dtaudevpar = BF0
    if time_option == 1
    # here the integration time is the proper time (simplest case)
        dtaudevpar = BF1
    elseif time_option == 2
    # here the integration time is an over increasing function of proper time (see example below,
    # where tau is given as a Taylor series wrt to proper time with small positive coefficients)
    # dtaudevpar must be ideally between 0 and 1
        dtaudevpar_tau0 = BF1-BF1/BigFloat(24327) # for test, June 25, 2024, here the integration time is approximatively running at the same rate than the proper time
        d2taudevpar2_tau0 = -BF1/BigFloat(34572) # # for test, June 25
        d3taudevpar3_tau0 = -BF1/BigFloat(45628) # # for test, June 25
        #tau = tau0+dtaudevpar_tau0*(evpar-evpar0)+d2taudevpar2_tau0*(evpar-evpar0)^2/BF2+d3taudevpar3_tau0*(evpar-evpar0)^3/BF6 # Taylor's series, this is an example of relationship, function must be monotonic, always increasing function of each other
        dtaudevpar = dtaudevpar_tau0+d2taudevpar2_tau0*(evpar-evpar0)+d3taudevpar3_tau0*(evpar-evpar0)^2/BF2  # first derivative of the previous line wrt integration time
        #d2taudevpar2 = d2taudevpar2_tau0+d3taudevpar3_tau0*(evpar-evpar0) # same, but for second derivative
        if dtaudevpar < BF0
            prt3("negative derivative of proper time wrt integration time (time_option = 2)")
            throw(error)
        end
        if dtaudevpar > BF1
            prt3("warning: derivative of proper time wrt integration time (time_option = 2) larger than 1")
        end
    elseif  time_option == 3
    # here the integration time is the coordinate time. In this case dtaudevpar and d2taudevpar2 are computed by imposing
    # that the partial second derivative of the coordinate time wrt to the coordinate time is zero
    # to be implemented, or probably rewritten, as it will be pretty inefficient here
        dtaudevpar = BF1/y[2]
        if dtaudevpar < BF0
            prt3("negative derivative of proper time wrt integration time (time_option = 3)")
            throw(error)
        end
        if dtaudevpar > BF1
            prt3("warning: derivative of proper time wrt integration time (time_option = 3) larger than 1")
        end
    elseif time_option == 4
        dtcdevpar_tc0 = BF1-BF1/BigFloat(24627) # for test, June 25, 2024, here the integration time is approximatively running at the same rate than the proper time
        d2tcdevpar2_tc0 = -BF1/BigFloat(34272) # # for test, June 25
        d3tcdevpar3_tc0 = -BF1/BigFloat(45728) # # for test, June 25
        #example of relationship, function must be monotonic, always increasing function of each other
        dtcdevpar = dtcdevpar_tc0+d2tcdevpar2_tc0*(evpar-evpar0)+d3tcdevpar3_tc0*(evpar-evpar0)^2/BF2  # first derivative of the previous line wrt integration time
        dtaudevpar = BF1/y[2]*dtcdevpar
        if dtaudevpar < BF0
            prt3("negative derivative of proper time wrt integration time (time_option = 4)")
            throw(error)
        end
        if dtaudevpar > BF1
            prt3("warning: derivative of proper time wrt integration time (time_option = 4) larger than 1")
        end
    else
        prt3("non-existing time scale for time_option")
        throw(error)
    end
    return dtaudevpar
end
function der1_tau_internal(y::Vector{BigFloat}, Given_Metric::Function)
# computes internally the numerical first derivative of proper time tau wrt to the integration time evpar
# IF and only IF the velocity is expressed in chart coordinates (not Cholesky of FW tetrad)
# this should be 1 numerically iff the integration time is the proper time
    (mt_s, ) = Given_Metric(y) # covariant and contravariant metric    
    dsdevpar2 = BF0 
    for p = 1:4
        for q = 1:4
            dsdevpar2 = dsdevpar2+mt_s[p, q]*y[2*p]*y[2*q]
        end
    end
    dtaudevpar = sqrt(dsdevpar2)/c_light # this is relative to the 4-velocity wrt the integration time evpar
    return dtaudevpar
end
function Numerical_parder_metric_y(order_parder_metric::Int, y::Vector{BigFloat}, Given_Metric::Function)
# partial derivatives of the metric, computed w.r.t order_parder_metric defined in main
    Permitted_order = Set([2, 4, 6, 8])
    if !(order_parder_metric in Permitted_order)
        throw(ErrorException)
    end
    parder_metric_y = fill(BF0, (4, 4, 4)) # first index = index of coordinate derivative
    # nsyst = length(y)
    cdc = fill(BF0, (4, 4)) # central difference coefficients
    index = div(order_parder_metric, 2)
    cdc[1, 1] = BF1/BF2				
    cdc[2, 1] = BF2/BF3 ; cdc[2, 2] = −BF1/BF12			
    cdc[3, 1] = BF3/BF4 ; cdc[3, 2] = −BF3/BF20 ; cdc[3, 3] = BF1/BF60		
    cdc[4, 1] = BF4/BF5 ; cdc[4, 2] = −BF1/BF5  ; cdc[4, 3] = BF4/BF105 ; cdc[4, 4] = −BF1/BF280
    for i = 1:4
        eps_CSN = max(epsmin_CSN, epsabs_CSN*abs(y[2*i-1])) # epsmin_CSN, epsabs_CSN defined in main, this must be checked to be sure that these values are OK for the case being studied
        ycall = deepcopy(y)
        for p = 1:index
            ycall[2*i-1] = y[2*i-1]+BigFloat(p)*eps_CSN
            (mt_p, ) = Given_Metric(ycall)
            for j = 1:4
                for k = 1:4   
                    parder_metric_y[i, j, k] = parder_metric_y[i, j, k]+cdc[index, p]*mt_p[j, k]/eps_CSN
                end
            end
            ycall[2*i-1] = y[2*i-1]-BigFloat(p)*eps_CSN
            (mt_m, ) = Given_Metric(ycall)
            for j = 1:4
                for k = 1:4
                    parder_metric_y[i, j, k] = parder_metric_y[i, j, k]-cdc[index, p]*mt_m[j, k]/eps_CSN
                end
            end
        end
    end
    return parder_metric_y
end
function Numerical_Christoffel_Symbols(order_parder_metric::Int, y::Vector{BigFloat}, Given_Metric::Function)
    # Build numerically Christoffel Symbols
    BigGamma = fill(BF0, (4, 4, 4)) # first indice = upper position of the Christoffel symbol, the two lower indices are symetrical
    parder_metric_y = Numerical_parder_metric_y(order_parder_metric, y, Given_Metric)
    (mt, mt_i) = Given_Metric(y)
    for i = 1:4
        for j = 1:4
            for k = 1:4
                s = BF0
                for p = 1:4
                s = s+mt_i[i, p]*(parder_metric_y[j, p, k]+parder_metric_y[k, j, p]-parder_metric_y[p, j, k])
                end
                BigGamma[i, j, k] = s/BF2
            end
        end
    end
    return BigGamma
end
function dot_product_at_y(vector1::Vector{BigFloat}, vector2::Vector{BigFloat}, Given_Metric::Function, y::Vector{BigFloat})
    # returns the dot product of two vectors wrt the metric at the chart point contained in y
    (mt_s, ) = Given_Metric(y) # covariant and contravariant metric
    dot_product = BF0
    for i = 1:4
        for j = 1:4
            dot_product = dot_product+mt_s[i, j]*vector1[i]*vector2[j]
        end
    end
    return dot_product
end
function norm2_at_y(vector::Vector{BigFloat}, Given_Metric::Function, y::Vector{BigFloat})
# to verify, if needed, the nature of a 4-vector. This is an invariant wrt to a change of coordinates
    norm2 = dot_product_at_y(vector, vector, Given_Metric, y)
    if norm2 > BF0
        prt3("time-like vector = ", norm2)
     elseif norm2 < BF0
        prt3("space-like vector = ", norm2)
     else
        prt3("null-vector = ", norm2)
     end
    return norm2
    end
function Vierbein(SetLinIndVect::Matrix{BigFloat}, Given_Metric::Function, y::Vector{BigFloat})
# SetLinIndVect :: Set of 4 *Linearly Independent* Vectors of strictly positive or strictly negative norms
# Ouput an orthonormal tetrad with respect to the Minkowski Cartesian metric (+c2, -1, -1, -1)
# Gram-Schmidt orthogonalization, then normalization (c, -1)
# if one or more of the provided vectors is a linear combination of the other vectors,
# the routine will fail with likely outputted NaN
# be aware of the confusion between the names of tetrad and vierbein, due to the definition of a vierbein
# in in this routine, two matrices are outputted, tetrad and tetrad_inv, that are in fact the vierbein and inverse vierbein transforms
# according to: chart_vector = tetrad*tetrad_vector and tetrad_vector = tetrad_inv*chart_vector
    (mt_s, mt_si) = Given_Metric(y) # covariant and contravariant metric
    tetrad = fill(BF0, (4, 4))
    # Gram-Schmidt construction
    for k = 1:4
        for i = 1:4
            tetrad[i, k] = SetLinIndVect[i, k]
        end
        for n = 2:k
            NUM = BF0
            DEN = BF0
            for p = 1:4
                for q = 1:4
                    NUM = NUM+mt_s[p, q]*SetLinIndVect[p, k]*tetrad[q, n-1]
                    DEN = DEN+mt_s[p, q]*tetrad[p, n-1]*tetrad[q, n-1]
                end
            end
            for i = 1:4
                tetrad[i, k] = tetrad[i, k]-NUM/DEN*tetrad[i, n-1]
            end
        end
    end 
    # orthonormalization
    for k = 1:4
        norm2 = BF0
        for i = 1:4
            for j = 1:4
                norm2 = norm2+mt_s[i, j]*tetrad[i, k]*tetrad[j, k]
            end
        end
        norm = sqrt(abs(norm2)) # abs(norm2), because normalization means +/- 1, not only 1
        for i = 1:4
            tetrad[i, k] = tetrad[i, k]/norm
        end
    end
    for i = 1:4
        tetrad[i, 1] = tetrad[i, 1]*c_light
    end
# verification of validity of tetrad, just a verification of the definition of a tetrad
    tetrad_verif = fill(BF0, (4, 4))
    for i = 1:4
        for j = 1:4
            tetrad_verif[j, i] = BF0
            for p = 1:4
                for q = 1:4
                    tetrad_verif[j, i] = tetrad_verif[j, i]+mt_s[p, q]*tetrad[p, i]*tetrad[q, j]
                end
            end
        end
    end
    dif2 = BF0
    (mt_m, ) = Minkowski_Metric_Cartesian(y)
    for i = 1:3
        for j = 1:4
            dif2 = dif2+(tetrad_verif[i, j]-mt_m[i, j])^2
        end
    end
    dif = sqrt(dif2)
    if dif > eps_tetrad
        prt3("possible instability in Vierbein routine = ", Float64(dif))
    end
    dif = verif_tetrad(tetrad, y, Given_Metric)
    tetradinv = inv(tetrad)
return (tetrad, tetradinv)
end
function verif_tetrad(tetrad::Matrix{BigFloat}, y::Vector{BigFloat}, Given_Metric::Function)
# verifies if a proposed tetrad is really a tetrad, up to a given accuracy (eps_tetrad)
    (mt_s, ) = Given_Metric(y) # covariant and contravariant metric
    (mt_m, ) = Minkowski_Metric_Cartesian(y) # covariant and contravariant metric
    # verification of validity of tetrad
    tetrad_v = fill(BF0, (4, 4))
    for i = 1:4
        for j = 1:4
            tetrad_v[j, i] = BF0
            for p = 1:4
                for q = 1:4
                    tetrad_v[j, i] = tetrad_v[j, i]+mt_s[p, q]*tetrad[p, i]*tetrad[q, j]
                end
            end
        end
    end
    dif2 = BF0
    for i = 1:3
        for j = 1:4
            dif2 = dif2+(tetrad_v[i, j]-mt_m[i, j])^2
        end
    end
    dift = sqrt(dif2)
    if dift > eps_tetrad
        prt3("possible inaccuracy in tetrad = ", Float64(dift))
    end
    tetradinv = inv(tetrad)
    tetrad_v = fill(BF0, (4, 4))
    for i = 1:4
        for j = 1:4
            tetrad_v[j, i] = BF0
            for p = 1:4
                for q = 1:4
                    tetrad_v[j, i] = tetrad_v[j, i]+mt_m[p, q]*tetradinv[p, i]*tetradinv[q, j]
                end
            end
        end
    end
    dif2 = BF0
    for i = 1:3
        for j = 1:4
            dif2 = dif2+(tetrad_v[i, j]-mt_s[i, j])^2
        end
    end
    difti = sqrt(dif2)
    if difti > eps_tetrad
        prt3("possible inaccuracy in tetradinv = ", Float64(difti))
    end
    return (dift, difti)
end
function FWtetrad_from_y(y::Vector{BigFloat})
    # returns the tetrad/vierbein contained in state vector y, as well as the inverse tetrad
    # according to: chart_vector = tetrad*tetrad_vector and tetrad_vector = tetrad_inv*chart_vector
    tetrad = fill(BF0, (4, 4)) 
    for i = 1:4 # tetrad vector
        for j = 1:4 # vector elements
            tetrad[j, i] = y[4*i+j+6]
        end
    end
    tetradinv = inv(tetrad)
    return (tetrad, tetradinv)
 end
function y_from_FWtetrad(y::Vector{BigFloat}, tetrad::Matrix{BigFloat})
# modify the tetrad/vierbein part of y   
# according to: chart_vector = tetrad*tetrad_vector and tetrad_vector = tetrad_inv*chart_vector
    for i = 1:4 # tetrad vector
        for j = 1:4 # vector elements
            y[4*i+j+6] = tetrad[j, i]
        end
    end
    return y
end
function Cholesky_tetrad_der(y::Vector{BigFloat}, Given_Metric::Function)
    # returns: ((tetrad, tetradder = directional derivative of the tetrad along the path, tetradparder = partial derivatives of the tetrad along the path)
    # and same for the inverse of the tetrad
    # Cholesky-like tetrad, with only 10-free parameters, instead of the full 16 free parameters
    # can always be defined for any smooth metric
    tetradinv = fill(BF0, (4, 4))
    a = fill(BF0, (4, 4))
    M = fill(BF0, (4, 4))
    M[1, 1] = c_light^2 ; M[2, 2] = -BF1 ; M[3, 3] = -BF1 ; M[4, 4] = -BF1
    (mt_s, mt_si) = Given_Metric(y) # covariant and contravariant metric
    c2 = c_light^2
    a[1, 1] = sqrt(mt_s[1, 1]/c2) ; a[1, 2] = mt_s[1, 2]/a[1, 1]/c2 ; a[1, 3] = mt_s[1, 3]/a[1, 1]/c2 ; a[1, 4] = mt_s[1, 4]/a[1, 1]/c2
    a[2, 2] = sqrt(a[1, 2]^2*c2-mt_s[2, 2]) ; a[2, 3] = (a[1, 2]*a[1, 3]*c2-mt_s[2, 3])/a[2, 2] ; a[2, 4] = (a[1, 2]*a[1, 4]*c2-mt_s[2, 4])/a[2, 2]
    a[3, 3] = sqrt(a[1, 3]^2*c2-a[2, 3]^2-mt_s[3, 3]) ; a[3, 4] = (a[1, 3]*a[1, 4]*c2-a[2, 3]*a[2, 4]-mt_s[3, 4])/a[3, 3]
    a[4, 4] = sqrt(a[1, 4]^2*c2-a[2, 4]^2-a[3, 4]^2-mt_s[4, 4])
    tetradinv[1, 1] = a[1, 1] ; tetradinv[1, 2] = a[1, 2] ; tetradinv[1, 3] = a[1, 3] ; tetradinv[1, 4] = a[1, 4]
    tetradinv[2, 2] = a[2, 2] ; tetradinv[2, 3] = a[2, 3] ; tetradinv[2, 4] = a[2, 4]
    tetradinv[3, 3] = a[3, 3] ; tetradinv[3, 4] = a[3, 4]
    tetradinv[4, 4] = a[4, 4]
    tetrad = inv(tetradinv)
    b = fill(BF0, (10, 10))
    b[1, 1] =BF2*c2*a[1, 1]
    b[5, 2] =BF2*c2*a[1, 2]   ; b[5, 5] =-BF2*a[2, 2]
    b[8, 3] =BF2*c2*a[1, 3]   ; b[8, 6] =-BF2*a[2, 3] ; b[8, 8] =-BF2*a[3, 3]
    b[10, 4] = BF2*c2*a[1, 4]   ; b[10, 7] = -BF2*a[2, 4] ; b[10, 9] = -BF2*a[3, 4] ; b[10, 10] = -BF2*a[4, 4]
    b[2, 1] =c2*a[1, 2]       ; b[2, 2] =c2*a[1, 1]
    b[3, 1] =c2*a[1, 3]       ; b[3, 3] =c2*a[1, 1]
    b[4, 1] =c2*a[1, 4]       ; b[4, 4] =c2*a[1, 1]
    b[6, 2] =c2*a[1, 3]       ; b[6, 3] =c2*a[1, 2]   ; b[6, 5] = -a[2, 3]      ;  b[6, 6] = -a[2, 2]
    b[7, 2] =c2*a[1, 4]       ; b[7, 4] =c2*a[1, 2]   ; b[7, 5] = -a[2, 4]      ;  b[7, 7] = -a[2, 2]
    b[9, 3] =c2*a[1, 4]       ; b[9, 4] =c2*a[1, 3]   ; b[9, 6] = -a[2, 4]      ;  b[9, 7] = -a[2, 3] ; b[9, 8] = -a[3, 4] ; b[9, 9] = -a[3, 3]
    b[10, 4] = BF2*c2*a[1, 4]   ; b[10, 7] = -BF2*a[2, 4] ; b[10, 9] = -BF2*a[3, 4] ;  b[10, 10] = -BF2*a[4, 4]
    parder_metric_y = Numerical_parder_metric_y(order_parder_metric, y, Given_Metric)
    s_member = fill(BF0, (10, 4))
    LINE = 0
    for k = 1:4
        for l = k:4
            LINE = LINE+1
            for r = 1:4
                s_member[LINE, r] = s_member[LINE, r]+parder_metric_y[r, k, l]
            end
        end
    end
    sol = inv(b)*s_member
    tetradparder = fill(BF0, (4, 4, 4))
    tetradinvparder = fill(BF0, (4, 4, 4))
    for r = 1:4
        tetradinvparder[r, 1, 1] = sol[1, r] ; tetradinvparder[r, 1, 2] = sol[2, r] ; tetradinvparder[r, 1, 3] = sol[3, r] ; tetradinvparder[r, 1, 4] = sol[4, r]
        tetradinvparder[r, 2, 2] = sol[5, r] ; tetradinvparder[r, 2, 3] = sol[6, r] ; tetradinvparder[r, 2, 4] = sol[7, r]
        tetradinvparder[r, 3, 3] = sol[8, r] ; tetradinvparder[r, 3, 4] = sol[9, r]
        tetradinvparder[r, 4, 4] = sol[10, r]
    end
    tetradinvder = fill(BF0, (4, 4))
    for i = 1:4
        for j = 1:4
            for r = 1:4
                tetradinvder[i, j] = tetradinvder[i, j]+tetradinvparder[r, i, j]*y[2*r]
            end
        end
    end
    tetradder = -tetrad*tetradinvder*tetrad
    intermediate = fill(BF0, (4, 4))
    for r = 1:4
        for i = 1:4
            for j = 1:4
                intermediate[i, j] = tetradinvparder[r, i, j]
            end
        end
        intermediate = -tetrad*intermediate*tetrad
        for i = 1:4
            for j = 1:4
                tetradparder[r, i, j] = intermediate[i, j]
            end
        end
    end
    # tetradparder and tetradinvparder can be returned if necessary
    return ((tetrad, tetradder), (tetradinv, tetradinvder))
end
function pieces_of_y(y::Vector{BigFloat})
    # pieces of y (or yp)
    # the velocities are either wrt proper time / integration time / global chart / local tetrad, that does not matter
    position = fill(BF0, 4)
    velocity = fill(BF0, 4)
    coordinate_time = fill(BF0, 2)
    FW_tetrad = fill(BF0, (4, 4))
    for i = 1:4
        position[i] = y[2*i-1]
        velocity[i] = y[2*i]
    end
    coordinate_time[1] = y[9]
    coordinate_time[2] = y[10]
    for i = 1:4 # tetrad vector
        for j = 1:4 # elements of each vector tetrad, each column vector = a tetrad vector
            FW_tetrad[j, i] = y[4*i+j+6]
        end
    end
    return (position, velocity, coordinate_time, FW_tetrad)
end
function y_from_pieces(position::Vector{BigFloat}, velocity::Vector{BigFloat}, coordinate_time::Vector{BigFloat}, FW_tetrad::Matrix{BigFloat})
    # y (or yp) from pieces of y (or yp)
    # nsyst should be = length(position)+lenght(velocity)+length(coordinate_time)+length(FW_tetrad)
    # the velocities are either wrt proper time / integration time / global chart / local tetrad, that does not matter
    y = fill(BF0, nsyst)
    for i = 1:4
        y[2*i-1] = position[i]
        y[2*i] = velocity[i]
    end
    y[9] = coordinate_time[1]
    y[10] = coordinate_time[2]
    for i = 1:4 # tetrad vector
        for j = 1:4 # elements of each vector tetrad, each column vector = a tetrad vector
            y[4*i+j+6] = FW_tetrad[j, i]
        end
    end
    return y
end
function velocity_from_y(y::Vector{BigFloat})
    # extract the velocity vector from y
    velocity = fill(BF0, 4)
    for p = 1:4
        velocity[p] = y[2*p]
    end
    return velocity
    end
function chart_vector_from_tetrad_vector(local_vector::Vector{BigFloat}, y::Vector{BigFloat})
    # from vector components wrt to tetrad frame to vector components wrt the local holonomic frame
    (tetrad, tetradinv) = FWtetrad_from_y(y)
    chart_vector = fill(BF0, 4)
    for p = 1:4
        for q = 1:4
            chart_vector[p] = chart_vector[p]+tetrad[p, q]*local_vector[q]
        end
    end
    return chart_vector
end
function tetrad_vector_from_chart_vector(chart_vector::Vector{BigFloat}, y::Vector{BigFloat})
    # from vector components wrt the local holonomic frame to vector components wrt to tetrad frame
    (tetrad, tetradinv) = FWtetrad_from_y(y)
    tetrad_vector = fill(BF0, 4)
    for p = 1:4
        for q = 1:4
            tetrad_vector[p] = tetrad_vector[p]+tetradinv[p, q]*chart_vector[q]
        end
    end
    return tetrad_vector
end
function ychol_from_ychart(y::Vector{BigFloat}, Given_Metric)
    # position given in chart coordinates,
    # velocity given in chart coordinates (chart tetrad to be totally correct), any time of integration
    (position_in_chart, velocity_in_chart, dummy1, dummy2) = pieces_of_y(y)
    ((tetradchol, tetradcholder), (tetradcholinv, tetradcholinvder)) = Cholesky_tetrad_der(y, Given_Metric)
    velocity_in_chol_tetrad = tetradcholinv*velocity_in_chart # from velocity in chart to local Cholesky tetrad
    ychol = y_from_pieces(position_in_chart, velocity_in_chol_tetrad, dummy1, dummy2)
    return ychol
end
function ychart_from_ychol(ychol::Vector{BigFloat}, Given_Metric)
    # position given in chart coordinates,
    # velocity given in Cholesky tetrad, any time of inegration
    (position_in_chart, velocity_in_chol_tetrad, dummy1, dummy2) = pieces_of_y(ychol)
    # only position (so in chart coordinates) is used inside Cholesky_tetrad_der, so we can input ychol in Cholesky_tetrad_der
    ((tetradchol, tetradcholder), (tetradcholinv, tetradcholinvder)) = Cholesky_tetrad_der(ychol, Given_Metric)
    velocity_in_chart = tetradchol*velocity_in_chol_tetrad
    chart_y = y_from_pieces(position_in_chart, velocity_in_chart, dummy1, dummy2)
    return chart_y
end
function ypchol_from_ypchart(y::Vector{BigFloat}, yp::Vector{BigFloat}, Given_Metric)
    # velocity in chart coordinates
    # acceleration in chart coordinates
    #(position_in_chart, velocity_in_chart0, dmmy1, dumm2) = pieces_of_y(y::Vector{BigFloat})
    (velocity_in_chart, acc_in_chart, dummy1, dummy2) = pieces_of_y(yp::Vector{BigFloat})
    # caution: tetradcholinvder is wrt the time usd for the time derivative
    ((dummy3, dummy4), (tetradcholinv, tetradcholinvder)) = Cholesky_tetrad_der(y, Given_Metric)
    acc_in_chol_tetrad = tetradcholinvder*velocity_in_chart+tetradcholinv*acc_in_chart
    ypchol = y_from_pieces(velocity_in_chart, acc_in_chol_tetrad, dummy1, dummy2)
    return ypchol
end
function ypchart_from_ypchol(ychol::Vector{BigFloat}, ypchol::Vector{BigFloat}, Given_Metric)
    #(position_in_chart, velocity_in_chart0, dmmy1, dumm2) = pieces_of_y(ychol::Vector{BigFloat})
    yt = ychart_from_ychol(ychol::Vector{BigFloat}, Given_Metric)
    (velocity_in_chart, acc_in_chol_tetrad, dummy1, dummy2) = pieces_of_y(ypchol::Vector{BigFloat})
    ((tetradchol, tetradcholder), (tetradcholinv, dummy3)) = Cholesky_tetrad_der(yt, Given_Metric)
    # caution: tetradcholder is wrt the time usd for the time derivative
    acc_in_chart = tetradcholder*tetradcholinv*velocity_in_chart+tetradchol*acc_in_chol_tetrad
    ypchart = y_from_pieces(velocity_in_chart, acc_in_chart, dummy1, dummy2)
    return ypchart
end
function Non_Grav_Grad(y::Vector{BigFloat})
    # non-gravitational gradient (a tensor) of the energy-momentum tensor, as only a function of 4-chart coordinates
    # only the part that is orthogonal to the worldline is seen by the spacecraft
    NGG = fill(BF0, 4) # contravariant tensor with respect to proper time
    # here an example for Parker Solar probe with actual values (radiation pressure from the Sun)
    NGG[1] = BF0
    # const mass_PSP = parse(BigFloat, "655.0") # mass in kg, PSP
    # const shield_area = parse(BigFloat, "4.0") # Sun shield area in m^2, PSP
    # const shiel_reflection_coef = parse(BigFloat, "1.8") # shield reflection coef, shield always facing the Sun at periapsis
    # radiation along the r polar chart coordinate here
    NGG[2] = solar_radiation_pressure*shield_area*shield_reflection_coef*(AU/y[3])^2/mass_PSP
    NGG[3] = BF0 # 
    NGG[4] = BF0 #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #NGG = fill(BF0, 4) # put zero for non-gravitional fordes, so no NG forces if this line is active
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return NGG
end
function Non_Grav_Bus(y::Vector{BigFloat})
    # non-gravitational forces measured on te spacecraft bus 
    # the first component IS, by definition, ALWAYS zero 
    NGB = fill(BF0, 4)
    NGB[1] = BF0 # DO NOT CHANGE IT, THIS COMPONENT IS ALWAYS ZERO BY ***DEFINITION*** IN THE S/C BUS
    NGB[2] = BF0 # as you want (model or measurement by S/C onboard accelerometers)
    NGB[3] = BF0 # as you want (idem)
    NGB[4] = BF0 # as you want (idem)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NGB = fill(BF0, 4) # put zero for the tests (Oct 16, 2024), so no NG forces here if this line is active
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return NGB
end
function Eqs_Motion_chol(evpar::BigFloat, ychol::Vector{BigFloat}, Given_Metric::Function, Non_Grav_Grad::Function, Non_Grav_Bus::Function)
    # the spherical coordinates integrated here are non-observable quantities
    # K contravariant gradient of energy-momentum tensor
    # the spherical coordinates integrated here are non-observable quantities
    # contravariant coordinates and velocities in spatial polar coordinates here
    # t = y[1] = coordinate time (clock at infinity) ; r = y[3] = pseudo-radius
    # phi = y[5] = pseudo-longitude ; theta = y[7] = pseudo-colatitude,
    # tau is proper time
    nsyst = length(ychol)
    ypchol = fill(BF0, nsyst)
    yp = fill(BF0, nsyst)
    # going back to y from ychol for the computations in Eqs_Motion_Chol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ychart = ychart_from_ychol(ychol, Given_Metric) # so in the integration time
    chart_four_velocity = velocity_from_y(ychart) # here, it is 4-velocity, i.e. wrt to proper time
    chol_velocity = velocity_from_y(ychol) # here it is velocity from ychol, wrt to the time of integration
    ((tetradchol, tetradcholder), (tetradcholinv, tetradcholinvder)) = Cholesky_tetrad_der(ychart, Given_Metric) # so here wrt integration time
    dtaudevpar_external = evpar2tau_external(time_option, evpar, ychart, Given_Metric) # this give from an external a priori the values of tau and its derivative wrt evpar
    ##-------------------------------------------------------------------------------------------------
    dtaudevpar_eqm = dtaudevpar_external
    ##-------------------------------------------------------------------------------------------------
    (mt_s, mt_si) = Given_Metric(ychart) # covariant and contravariant metric
    NGG = Non_Grav_Grad(ychart) # covariant gradient (provided) of the energy-momentum tensor
    NGB = Non_Grav_Bus(ychart) # Non-Gravitational forces measured wrt the spacecraft
    # order_parder_metric is defined in main
    BigGamma = Numerical_Christoffel_Symbols(order_parder_metric, ychart, Given_Metric)
    F = fill(BF0, 4) # contravariant 4-force, orthogonal to the wordline
    for i = 1:4
        for j = 1:4
            F[i] = F[i]+NGG[j]*(mt_si[j, i]-chart_four_velocity[j]*chart_four_velocity[i]/c_light^2)
        end
    end
    ypBG = fill(BF0, 4)
    yp_tot = fill(BF0, 4)
    yp_pt = fill(BF0, nsyst)
    for i = 1:4
        for j = 1:4
            for k = 1:4
                ypBG[i] = ypBG[i]+BigGamma[i, j, k]*chart_four_velocity[j]*chart_four_velocity[k] # Big Gamma gravity acceleration
            end
        end
    end
    for i = 1:4
        yp_tot[i] = ypBG[i]-F[i] # sum of Big Gamma gravity acceleration + non-gravitational forces (with a minus sign)
        yp_pt[2*i] = yp_tot[i]
    end
    ypchol_vel = -dtaudevpar_eqm*tetradcholinv*(tetradcholder*chol_velocity+yp_tot) # as ypBG and F are evaluated wrt proper time, tetradcholder wet to integration time (ychart)
    for i = 1:4
        ypchol[2*i-1] = dtaudevpar_eqm*ychart[2*i] # careful, this is because we need to have the chart coordinates with respect to integration time
        ypchol[2*i] = ypchol_vel[i]
    end
    ypchol[9] = dtaudevpar_eqm # normally, dtaudevpar_external, so will give proper time by integration
    ypchol[10] = BF1 # integration time (elapsed)
# ------------------------------------------------computation of the transported tetrad---------------------------------------------------------------------------
# parallel transport wrt integration time of the spatial part of the tetrad
# caution: this will only work as a Fermi-walker transport if and only if the first vector of the tetrad is the 4-velocity
# otherwise the transported tetrad will rotate / wrt distant stars
    for p = 1:4
        for i = 1:4
            yp[4*p+i+6] = BF0
            for j = 1:4
                for k = 1:4
                    yp[4*p+i+6] = yp[4*p+i+6]-BigGamma[i, j, k]*ychart[4*p+j+6]*ychart[2*k]
                end
            end
        end
    end
# Khl = correction to the parallel transport to get the Fermi-Walker transport
    Khh = fill(BF0, (4, 4))
    for p = 1:4
        for q = 1:4
            Khh[p, q] = (F[p]*ychart[2*q]-F[q]*ychart[2*p])/c_light^2
        end
    end
    Khl = fill(BF0, (4, 4))
    for p = 1:4
        for q = 1:4
            for l = 1:4
                Khl[p, q] = Khl[p, q]+mt_s[l, q]*Khh[p, l]
            end
        end
    end
# application of the Khl correction to get the Fermi-Walker transport
    for p = 1:4
        for i = 1:4
            for l = 1:4
                yp[4*p+i+6] = yp[4*p+i+6]+Khl[i, l]*ychart[4*p+l+6]
            end
        end
    end
    for i = 11:26
        ypchol[i] = dtaudevpar_eqm*yp[i]
    end
    return ypchol
end

# ================================================================================================
# 🌌 Simplified relativistic equation of motion — no tetrad, no Cholesky
# ================================================================================================
function Eqs_Motion_chol_simple(evpar::BigFloat,
                                y::Vector{BigFloat},
                                Given_Metric::Function,
                                Non_Grav_Grad::Function,
                                Non_Grav_Bus::Function)

    dydσ = similar(y)

    # Metric (g, ginv) and Christoffel symbols Γ^μ_{αβ}
    (g, ginv) = Given_Metric(y)
    Γ = Numerical_Christoffel_Symbols(order_parder_metric, y, Given_Metric)

    # --- 4-position derivatives wrt integration parameter σ ---
    # Ici on suppose que y[2k] = dx^μ/dσ (cohérent avec l’EDO simple)
    dydσ[1] = y[2]   # dt/dσ
    dydσ[3] = y[4]   # dr/dσ
    dydσ[5] = y[6]   # dφ/dσ
    dydσ[7] = y[8]   # dθ/dσ

    # --- 4-velocity derivatives (geodesic): d²x^μ/dσ² = - Γ^μ_{αβ} (dx^α/dσ)(dx^β/dσ)
    u = (BigFloat[y[2], y[4], y[6], y[8]])
    for μ in 1:4
        acc = BF0
        for α in 1:4, β in 1:4
            acc -= Γ[μ, α, β] * u[α] * u[β]
        end
        dydσ[2*μ] = acc
    end

    # --- time bookkeeping slots ---
    # dτ/dσ imposé par le choix de l’échelle de temps (time_option)
    dydσ[9]  = evpar2tau_external(time_option, evpar, y, Given_Metric)  # dτ/dσ
    dydσ[10] = BF1                                                      # dσ/dσ = 1

    # --- unused tetrad slots (11..26) ---
    for i in 11:length(y)
        dydσ[i] = BF0
    end

    return dydσ
end

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
function rki(t0::BigFloat, t1::BigFloat, y0::Vector{BigFloat}, Given_Metric::Function, Eqs_Motion::Function, Non_Grav_Grad::Function, Non_Grav_Bus::Function)
    # apply to Eqs_Motion or Eqs_Motion_chol indifferently
    # Gauss 5-stages implicit symplectic integrator
    # Butcher's 1962 paper
        (stages, ) = size(a_B) ; nsyst = length(y0)
        # the k_B's are initialized to zero, pretty rough, but it works,
        # the radius of convergence seems to be large, to be asserted (future work)
        k_B = fill(BF0, (stages, nsyst))
        h_B = t1-t0
        # the implicit method is implemented by a pretty rough fixed point method,
        # with an a priori fixed number of loops (itmax_rki) 
        # itmax_rki = 10 for Parker Solar Probe, defined in main
        # to be implemented: better initialization, convergence criteria following
        # the Hairer et al. 2002 book, p 275, future work
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
        for it = 1:itmax_rki
            for i = 1:stages
                s_B = fill(BF0, (nsyst, 1))
                for j = 1:stages
                    for n = 1:nsyst
                        s_B[n] = s_B[n]+a_B[i, j]*k_B[j, n]
                    end
                end
                yf = Eqs_Motion(t0+h_B*c_B[i, 1], vec(y0+h_B*s_B), Given_Metric, Non_Grav_Grad, Non_Grav_Bus)
                for n = 1:nsyst
                    k_B[i, n] = yf[n]
                end
            end
        end 
        dy = (b_B*k_B)'
        y1 = vec(y0+h_B*dy)
        return y1
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
function main(evpar0::BigFloat, ychol0::Vector{BigFloat}, Given_Metric::Function)
    nsyst = length(ychol0)
    y0 = ychart_from_ychol(ychol0, Given_Metric)
    dtaudevpar_internal = BF0 # do not remove, needed for scope of variable 
    y1 = fill(BF0, nsyst) # do not remove, otherwise this parameter is only defined inside the next loop !!!!
    dtaudevpar = BF0 # do not remove, otherwise this parameter is only defined inside the next loop !!!!
    ev_par0 = evpar0 # just for integration # start of time step
    prt3("entering main loop at = ", now())
        open("SCephemeris_$(timetag).txt", "w") do f
            for i = 1:N_integration_steps
                for j = 1:nsyst
                    write(f, "$(j) "*"$(y0[j]) "*"\n")
                end     
                ev_par1 = ev_par0+evpar_step # end of time step
                ychol1 = rki(ev_par0, ev_par1, ychol0, Given_Metric, Eqs_Motion_chol, Non_Grav_Grad, Non_Grav_Bus) # integration in integration time, quantities given at time = evpar0
                y1 = ychart_from_ychol(ychol1, Given_Metric)
                dtaudevpar_internal = der1_tau_internal(y1, Given_Metric)
                ev_par0 = ev_par1 ; y0 = deepcopy(y1) ; ychol0 = deepcopy(ychol1)    
            end # i loop  
            for j = 1:nsyst
                write(f, "$(j) "*"$(y0[j]) "*"\n")
            end
        end # open
        prt3("exiting main loop at = ", now())
     return ev_par0, ychol0, dtaudevpar_internal # tuple of values returned by main()
end
