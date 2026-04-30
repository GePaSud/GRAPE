module GRAPE_P
# GRAPE = General Relativity Accelerometer-based Propagation Environment
# coding in Julia by Jean-Pierre Barriot, April 12, 2026
# subject to Apache 2.0 license, GitHub.
# submitted to software-X journal (R1)
# Windows 11 version, Julia 1.12.1
#------------------------------------------------------------------------------------------------------------------------------------------------
# integration of the motion of a spacecraft (S/C) in a general relativity framework with associated attitude by Fermi-Walker transport,
# in presence (or not) of non-gravitational forces
# non-gravitational forces are either coded wrt to chart coordinates, as the gradient of the energy-momentum tensor, or as locally modeled or measured 
# accelerations in the S/C bus frame (so the name Accelerometer-based)
#------------------------------------------------------------------------------------------------------------------------------------------------
# symplectic integrator with conservation of the norm of the 4-velocity
# Gauss 5-stages symplectic integrator (Butcher, Math. Comput. 18(85), 50-64, 1964, ISSN 0025-5718)
# symplecticity is enforced by integrating the equations of motion of the S/C wrt to a specific tetrad along the worldline of the S/C
# named Cholesky tetrad in the Julia code
# see O'Leary and Barriot, 2021, Celestial Mech. Dyn. Astro. https://doi.org/10.1007/s10569-021-10051-7, for the mathematics
# see J.-P. Barriot, J. O'Leary and J. Yan, Beyond Newtonian Orbitography for Geodesy, Astronomy and Planetary Sciences:
# the GRAPE Project, Journal of Physics: Conference Series 3109 (2025) 012057, doi:10.1088/1742-6596/3109/1/012057
# for specifically the algorithm of this code
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
# as the connections are computed numerically,
# with the restriction that the first coordinate of integration (and the metric) must be the time coordinate.
# The integration time can be either proper time or coordinate time, or a function of them (see routine evpar2tau_external)
#------------------------------------------------------------------------------------------------------------------------------------------------
# all floating points quantities are coded with **extended precision** using the BigFloat format, do **NOT** suppress this feature,
# it is absolutely needed for a proper modeling of the relativistic effects, see keyword
# setprecision(BigFloat,xxx) # xxx=53 for 64 bits IEEE754, 113 for 128 bits, 237 for 256 bits IEEE754,
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
# tests were done with Julia version 1.12.1 
#----------------------------------------------------------- code starts here------------------------------------------------------------------------------
# ------------------------------------------------- first some comfort routines ---------------------------------------------------------------------------
function Call_System(command::String)
# to call the underlying Windows shell # not difficult to modify on Linux/Mac computers
# windows OS only, most of the basic commands should be "cmd /c ....", /c is to exit shell
# command="cmd /c cd $working_directory"
# system_message=Call_System(command)
# println("system_message= ",system_message)
        if !isascii(command)
           error("input of Call_System must be an all-ascii string")
        end
        split_string=split(command)
        return_string=run(`$split_string`)
    return return_string
end
function prt2(x...)
# Python-like print, with spaces printed between values for ease of use
    N=length(x)
    for i=1:N-1
        print(x[i]," ")
    end
    println(x[N])
    #write(runlogfile,timetag*"\n")
end
function prt3(x...)
# combines prt3 and prtstr to save in the runlogfile
    N=length(x)
    str=""
    for i=1:N-1
        print(x[i]," ")
        str=str*string(x[i])*" "
    end
    println(x[N])
    str=str*string(x[N])*"\n"
    write(runlogfile,str)
    return
end
function dumpcode(timetag::String)
# copy the running code into a file named codefile_timetag.txt
# global script is the name of the umbrella script running the GRAPE script (for example, an IDE debugger)
# local script is the name of the GRAPE script itself
# global script and local script names should be the same if julia is run from a terminal window
global_script=PROGRAM_FILE
prt3("global script=",global_script)
local_script=@__FILE__
prt3("local script=",local_script)
prt3("-------------dump of code being executed----------")
open(local_script,"r") do code # reading the code
codefile=open("codefile_$(timetag).txt","w")
    while !eof(code)
        line=readline(code)
        line=rstrip(lstrip(line))
        write(codefile,line*"\n")
    end
    close(codefile)
end
prt3("--------------end of dump of the code-------------")
end
function Give_Function_Name(Function_Name::Function)
    return Function_Name
end
# ------------------------------------------------- start of physics routines --------------------------------------------------------
function Minkowski_Metric_Cartesian(y::Vector{BigFloat})
# in Cartesian local coordinates, NOT in polar form
# y must be provided but is not used, to avoid a possible compilation error (empty entry list)
# Minkowski covariant metric tensor (c^2,-1,-1,-1).
# Note that c is given explicitly, and appears in the Metric, not in the coordinate definition
        mt_m=fill(BF0,(4,4))
        mt_m[1,1]=c_light^2
        mt_m[2,2]=-BF1 ; mt_m[3,3]=-BF1 ; mt_m[4,4]=-BF1
        # Minkowski contravariant metric tensor
        mt_mi=fill(BF0,(4,4))
        mt_mi[1,1]=BF1/mt_m[1,1] ; mt_mi[2,2]=BF1/mt_m[2,2]
        mt_mi[3,3]=BF1/mt_m[3,3] ; mt_mi[4,4]=BF1/mt_m[4,4]
        return (mt_m,mt_mi)
end
function Minkowski_Metric_Polar(y::Vector{BigFloat})
# for tests, polar coordinates
# t=y[1]=coordinate time (clock at infinity) ; r=y[3]=pseudo-radius
# phi=y[5]=pseudo-longitude ; theta=y[7]=pseudo-colatitude,
        r=y[3] ; theta=y[7]
        sn=sin(theta)
        # Minkowski covariant metric tensor
        mt_s=fill(BF0,(4,4))
        mt_s[1,1]=+c_light^2 ; mt_s[2,2]=-BF1          
        mt_s[3,3]=-(r*sn)^2    ; mt_s[4,4]=-r^2
        # Minkowski contravariant metric tensor
        mt_si=fill(BF0,(4,4))
        mt_si[1,1]=BF1/mt_s[1,1] ; mt_si[2,2]=BF1/mt_s[2,2]
        mt_si[3,3]=BF1/mt_s[3,3] ; mt_si[4,4]=BF1/mt_s[4,4]
        return (mt_s,mt_si)
end
function Schwarzschild_Metric_Polar(y::Vector{BigFloat})
# for tests, polar coordinates
# Schwarszchild metric, reduces to previous Minkowski_Metric_Cartesian if mass of star is zero (value rs=0)
# t=y[1]=coordinate time (clock at infinity) ; r=y[3]=pseudo-radius
# phi=y[5]=pseudo-longitude ; theta=y[7]=pseudo-colatitude,
    r=y[3] ; theta=y[7]
    sn=sin(theta)
    w=BF1-rs/r # rs = Schwarzschild radius
    # Schwarszchild covariant metric tensor
    mt_s=fill(BF0,(4,4))
    mt_s[1,1]=+c_light^2*w ; mt_s[2,2]=-BF1/w          
    mt_s[3,3]=-(r*sn)^2    ; mt_s[4,4]=-r^2
    # Schwarszchild contravariant metric tensor
    mt_si=fill(BF0,(4,4))
    mt_si[1,1]=BF1/mt_s[1,1] ; mt_si[2,2]=BF1/mt_s[2,2]
    mt_si[3,3]=BF1/mt_s[3,3] ; mt_si[4,4]=BF1/mt_s[4,4]
    return (mt_s,mt_si)
end
function Schwarzschild_BG_Polar(y::Vector{BigFloat})
# Christoffel symbols for Schwarszchild metric in polar coordinates, reduces to previous Minkowski_Metric_Cartesian if mass of star is zero (value rs=0)
# t=y[1]=coordinate time (clock at infinity) ; r=y[3]=pseudo-radius
# phi=y[5]=pseudo-longitude ; theta=y[7]=pseudo-colatitude,
# rs = Schwarzschild radius
    BG=fill(BF0,(4,4,4))
# indexes t ==1 ; r==2 ; 3 == longitude phi ; 4 == colatitude theta
    r=y[3] ; theta=y[7]
    BG[1,1,2]=BG[1,2,1]=rs/BF2/r/(r-rs)
    BG[2,1,1]=c_light^2*rs*(r-rs)/BF2/r^3
    BG[2,2,2]=-rs/BF2/r/(r-rs)
    BG[2,3,3]=-(r-rs)*(sin(theta))^2
    BG[2,4,4]=-(r-rs)
    BG[3,2,3]=BG[3,3,2]=BF1/r
    BG[3,3,4]=BG[3,4,3]=cot(theta)
    BG[4,3,3]=-sin(theta)*cos(theta)
    BG[4,2,4]=BG[4,4,2]=BF1/r
    return BG
end
function Schwarzschild_tetrad(y::Vector{BigFloat})
    (mt_s,mt_si)=Schwarzschild_Metric_Polar(y)
    tetrad=fill(BF0,(4,4))
    tetradinv=fill(BF0,(4,4))
    r=y[3] ; theta=y[7]
    sn=sin(theta)
    w=BF1-rs/r # rs = Schwarzschild radius
    tetradinv[1,1]=sqrt(w) ; tetradinv[2,2]=BF1/sqrt(w )         
    tetradinv[3,3]=r*sn    ; tetradinv[4,4]=r
    tetrad=inv(tetradinv)
    return(tetrad,tetradinv)
end
function TDB_FB(TT::BigFloat)
# TDB (coordinate, or ephemeris time) as a function of TT (Terrestrial Time)
# for illustration, not used in the test case example of SoftwareX
# Kopeikin et al. book (DOI:10.1002/9783527634569), page 715,
# abridged expression from Fairhead and Bretagnon, 1990
TT=Float64(TT)
T=(TT-2451545.0)/36525.0 # TT (Terrestrial Time) in JD
TDB=TT+0.001657  *sin( 628.3076*T+6.2401)
         +0.000022  *sin( 575.3385*T+4.2970)
         +0.000014  *sin(1256.6152*T+6.1969)
         +0.000005  *sin( 606.9777*T+4.0212)
         +0.000005  *sin(  52.9691*T+0.4444)
         +0.000002  *sin(  21.3299*T+5.5431)
         +0.000010*T*sin( 628.3076*T+4.2490)
TDB=parse(BigFloat,string(TDB))
    return TDB
end
function verify_harmonic_coord(y::Vector{BigFloat},Given_Metric::Function)
# verify if the coordinates are harmonic wrt the metric tensor
# for illustration, not used in the test case example of SoftwareX
# BigGamma # first indice = upper position of the Christoffel symbol,
# the two lower indices are symmetrical for holonomic coordinates (not mandatory for this code)
   (mt_c,mt_ci)=Given_Metric(y)
   BigGamma=Numerical_Christoffel_Symbols(order_parder_metric,y,Given_Metric)
   verif_harm=fill(BF0,4)
    for i=1:4
        for j=1:4
            for k=1:4
                verif_harm[i]=verif_harm[i]+mt_ci[j,k]*BigGamma[i,j,k]
            end
        end
    end
    return verif_harm
end
#=
# ------illustration, for possible use with Newtonian_Metric_Polar, not used for the test case example of SoftwareX-----------------------------------
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
   BODY_GM=fill(BF0,10) # units: KM**3/SEC**2
   BODY_GM[1]=parse(BigFloat,"22031.780000") # Mercury
   BODY_GM[2]=parse(BigFloat,"324858.592000") # Venus
   BODY_GM[3]=parse(BigFloat,"403503.235502") # Earth / Moon
   BODY_GM[4]=parse(BigFloat,"42828.375214") # Mars + moons
   BODY_GM[5]=parse(BigFloat,"126712764.800000") # Jupiter + moons
   BODY_GM[6]=parse(BigFloat,"37940585.200000") # Saturn + moons
   BODY_GM[7]=parse(BigFloat,"5794548.600000") # Uranus + moons
   BODY_GM[8]=parse(BigFloat,"6836527.100580") # Neptune + moons
   BODY_GM[9]=parse(BigFloat,"977.000000") # Pluto / Charon
   BODY_GM[10]=parse(BigFloat,"132712440041.939400") # Sun
#--------------------------------------------------------------------------------------
return BODY_GM
end
=#
function Newtonian_Metric_Polar(y::Vector{BigFloat})
    # this is the standard approximation of the Schwarschild metric by Newtonian potential, valid up to 1/c_light^4 order, with no frame-dragging
    r=y[3] ; theta=y[7]
    sn=sin(theta)
    newtonian_potential=mu/r # with a "+" sign, can be replaced by any valid potential model of the Solar System,
    # for exmple derived from DE430 (see previous lines, commented routine init_de430())
    mt_s=fill(BF0,(4,4))
    mt_s[1,1]=(c_light^2-BF2*newtonian_potential) ; mt_s[2,2]=-(BF1+BF2/c_light^2*newtonian_potential)        
    mt_s[3,3]=-(r*sn)^2    ; mt_s[4,4]=-r^2
    # Newtonian contravariant metric tensor
    mt_si=fill(BF0,(4,4))
    mt_si[1,1]=BF1/mt_s[1,1] ; mt_si[2,2]=BF1/mt_s[2,2]
    mt_si[3,3]=BF1/mt_s[3,3] ; mt_si[4,4]=BF1/mt_s[4,4]
    return (mt_s,mt_si)
# more complicated weak-field metrics, including PPN parameters,
# can be constructed here (see for example the book of Brumberg, https://doi.org/10.1201/9780203756591)
# and using the SPICE kernels: see https://naif.jpl.nasa.gov/pub/naif/misc/toolkit_docs_N0067/C/req/naif_ids.html
end
#----------this is the routine currently implemented in this code for the test case of SoftwareX-------------------
function Kerr_Metric_Polar(y::Vector{BigFloat})
# Kerr metric, reduces to Schwarzschild_Metric_Polar if the star is not rotating (value as=0)
# Boyer–Lindquist polar coordinates (technically oblate spheroidal coordinates),
# see catalog of spacetimes, https://doi.org/10.48550/arXiv.0904.4184
# the "as" parameter is the scaled angular momentum, in units of length
# as=Jstar/Mstar/c_light
    r=y[3] ; theta=y[7]
    sn=sin(theta) ; cs=cos(theta)
    BigSigma=r^2+as^2*cs^2
    BigDelta=r^2-rs*r+as^2
    mt_k=fill(BF0,(4,4))
    mt_k[1,1]=+c_light^2*(BF1-rs*r/BigSigma)
    mt_k[1,3]=mt_k[3,1]=rs*as*r*sn^2*c_light/BigSigma # multiplied by 2 by symmetry
    mt_k[2,2]=-BigSigma/BigDelta
    mt_k[3,3]=-(r^2+as^2+rs*as^2*r*sn^2/BigSigma)*sn^2
    mt_k[4,4]=-BigSigma
    mt_ki=inv(mt_k)
    return (mt_k,mt_ki)
    end
#function Any_SpaceTime_Metric(y)
# you can provide here the metric you want, including metrics build on DE430 JPL ephemeris or others, provided that the first coordinate
# MUST be the coordinate time (TBD, ephemeris_time, i.e. a clock at spatial infinity). This is because coordinate time is the 
# the time that links all events in the Solar System on a common clock
# signature must be (+,-,-,-)
#end
#
function Minkowski_norm2(vector::Vector{BigFloat})
# returns the Minkowskian (Cartesian) squared norm of a 4-vector
    norm2=c_light^2*vector[1]^2-vector[2]^2-vector[3]^2-vector[4]^2
    return norm2
end
function evpar2tau_external(time_option::Integer,evpar::BigFloat,y::Vector{BigFloat},Given_Metric::Function)
# gives the proper time as a function imposed by an a priori relationship, EXTERNAL to the code
# for example, the proper time of the spacecraft wrt to an Earth clock
# here tau=sigma, sigma being noted evpar (EVolution_PARameter) as variable in the julia code
# nor y nor Given_Metric are used in this routine as on April 12, 2026, reserved for future use
# a practical example linking terrestrial time (TT) with coordinate time TDB (or ET) can be found in Kopeikin et al.,, eq. 9.6, page 715
# at DOI:10.1002/9783527634569
    dtaudevpar=BF0
    if time_option == 1
    # here the integration time is the proper time (simplest case)
        dtaudevpar=BF1
    elseif time_option == 2
    # to be implemented by end-user
    # here the integration time is an over increasing function of proper time (see example below,
    # where tau is given as a Taylor series wrt to proper time with small positive coefficients)
    # dtaudevpar must be ideally between 0 and 1
        dtaudevpar_tau0=BF1-BF1/BigFloat(24327) # for example, here the integration time is approximatively running at the same rate than the proper time
        d2taudevpar2_tau0=-BF1/BigFloat(34572)  # # for example
        d3taudevpar3_tau0=-BF1/BigFloat(45628)  # # for example
        #tau=tau0+dtaudevpar_tau0*(evpar-evpar0)+d2taudevpar2_tau0*(evpar-evpar0)^2/BF2+d3taudevpar3_tau0*(evpar-evpar0)^3/BF6 # Taylor's series, this is an example of relationship, function must be monotonic, always increasing function of each other
        dtaudevpar=dtaudevpar_tau0+d2taudevpar2_tau0*(evpar-evpar0)+d3taudevpar3_tau0*(evpar-evpar0)^2/BF2  # first derivative of the previous line wrt integration time
        #d2taudevpar2=d2taudevpar2_tau0+d3taudevpar3_tau0*(evpar-evpar0) # same, but for second derivative
        if dtaudevpar < BF0
            prt3("negative derivative of proper time wrt integration time (time_option=2)")
            throw(error)
        end
        if dtaudevpar > BF1
            prt3("warning: derivative of proper time wrt integration time (time_option=2) larger than 1")
        end
    elseif  time_option == 3
    # here the integration time is the coordinate time. In this case dtaudevpar and d2taudevpar2 are computed by imposing
    # that the partial second derivative of the coordinate time wrt to the coordinate time is zero
        dtaudevpar=BF1/y[2]
        if dtaudevpar < BF0
            prt3("negative derivative of proper time wrt integration time (time_option=3)")
            throw(error)
        end
        if dtaudevpar > BF1
            prt3("warning: derivative of proper time wrt integration time (time_option=3) larger than 1")
        end
    elseif time_option == 4
        # to be implemented by end-user
        dtcdevpar_tc0=BF1-BF1/BigFloat(24627) # for example, here the integration time is approximatively running at the same rate than the proper time
        d2tcdevpar2_tc0=-BF1/BigFloat(34272)  # # for example
        d3tcdevpar3_tc0=-BF1/BigFloat(45728)  # # for example
        #example of relationship, function must be monotonic, always increasing function of each other
        dtcdevpar=dtcdevpar_tc0+d2tcdevpar2_tc0*(evpar-evpar0)+d3tcdevpar3_tc0*(evpar-evpar0)^2/BF2  # first derivative of the previous line wrt integration time
        dtaudevpar=BF1/y[2]*dtcdevpar
        if dtaudevpar < BF0
            prt3("negative derivative of proper time wrt integration time (time_option=4)")
            throw(error)
        end
        if dtaudevpar > BF1
            prt3("warning: derivative of proper time wrt integration time (time_option=4) larger than 1")
        end
    else
        prt3("non-existing time scale for time_option")
        throw(error)
    end
    return dtaudevpar
end
function der1_tau_internal(y::Vector{BigFloat},Given_Metric::Function)
# computes internally the numerical first derivative of proper time tau wrt the integration time evpar
# IF and only IF the velocity is expressed in chart coordinates (not Cholesky of FW tetrad)
# this should be 1 numerically if and only if the integration time is the proper time
    (mt_s,)=Given_Metric(y) # covariant and contravariant metric    
    dsdevpar2=BF0 
    for p=1:4
        for q=1:4
            dsdevpar2=dsdevpar2+mt_s[p,q]*y[2*p]*y[2*q]
        end
    end
    dtaudevpar=sqrt(dsdevpar2)/c_light # this is relative to the 4-velocity wrt the integration time evpar
    return dtaudevpar
end
function Numerical_parder_metric_y(order_parder_metric::Int,y::Vector{BigFloat},Given_Metric::Function)
# partial derivatives of the metric, computed w.r.t order_parder_metric defined in main
    Permitted_order=Set([2,4,6,8])
    if !(order_parder_metric in Permitted_order)
        throw(ErrorException)
    end
    parder_metric_y=fill(BF0,(4,4,4)) # first index = index of coordinate derivative
    # nsyst=length(y)
    cdc=fill(BF0,(4,4)) # central difference coefficients
    index=div(order_parder_metric,2)
    cdc[1,1]=BF1/BF2				
    cdc[2,1]=BF2/BF3 ; cdc[2,2]=−BF1/BF12			
    cdc[3,1]=BF3/BF4 ; cdc[3,2]=−BF3/BF20 ; cdc[3,3]=BF1/BF60		
    cdc[4,1]=BF4/BF5 ; cdc[4,2]=−BF1/BF5  ; cdc[4,3]=BF4/BF105 ; cdc[4,4]=−BF1/BF280
    for i=1:4
        eps_CSN=max(epsmin_CSN,epsabs_CSN*abs(y[2*i-1])) # epsmin_CSN,epsabs_CSN defined in main, this must be checked to be sure that these values are OK for the case being studied
        ycall=deepcopy(y)
        for p=1:index
            ycall[2*i-1]=y[2*i-1]+BigFloat(p)*eps_CSN
            (mt_p,)=Given_Metric(ycall)
            for j=1:4
                for k=1:4   
                    parder_metric_y[i,j,k]=parder_metric_y[i,j,k]+cdc[index,p]*mt_p[j,k]/eps_CSN
                end
            end
            ycall[2*i-1]=y[2*i-1]-BigFloat(p)*eps_CSN
            (mt_m,)=Given_Metric(ycall)
            for j=1:4
                for k=1:4
                    parder_metric_y[i,j,k]=parder_metric_y[i,j,k]-cdc[index,p]*mt_m[j,k]/eps_CSN
                end
            end
        end
    end
    return parder_metric_y
end
function Numerical_Christoffel_Symbols(order_parder_metric::Int,y::Vector{BigFloat},Given_Metric::Function)
    # Build numerically Christoffel Symbols
    BigGamma=fill(BF0,(4,4,4)) # first indice = upper position of the connection symbol,
    # the two lower indices are symmetrical (connection = Christoffel symbol in the case of holonomic coordinates)
    # in this case, the following computation can be shortened by taking into account this symmetry (not done here)
    parder_metric_y=Numerical_parder_metric_y(order_parder_metric,y,Given_Metric)
    (mt,mt_i)=Given_Metric(y)
    for i=1:4
        for j=1:4
            for k=1:4
                s=BF0
                for p=1:4
                s=s+mt_i[i,p]*(parder_metric_y[j,p,k]+parder_metric_y[k,j,p]-parder_metric_y[p,j,k])
                end
                BigGamma[i,j,k]=s/BF2
            end
        end
    end
    return BigGamma
end
function dot_product_at_y(vector1::Vector{BigFloat},vector2::Vector{BigFloat},Given_Metric::Function,y::Vector{BigFloat})
    # returns the dot product of two vectors wrt the metric at the chart point contained in y
    (mt_s,)=Given_Metric(y) # covariant and contravariant metric
    dot_product=BF0
    for i=1:4
        for j=1:4
            dot_product=dot_product+mt_s[i,j]*vector1[i]*vector2[j]
        end
    end
    return dot_product
end
function norm2_at_y(vector::Vector{BigFloat},Given_Metric::Function,y::Vector{BigFloat})
# to verify, if needed, the nature of a 4-vector. This is an invariant wrt a change of coordinates
    norm2=dot_product_at_y(vector,vector,Given_Metric,y)
    if norm2 > BF0
        prt3("time-like vector=",norm2)
     elseif norm2 < BF0
        prt3("space-like vector=",norm2)
     else
        prt3("null-vector=",norm2)
     end
    return norm2
    end
function Vierbein(SetLinIndVect::Matrix{BigFloat},Given_Metric::Function,y::Vector{BigFloat})
# SetLinIndVect :: Set of 4 *Linearly Independent* Vectors of strictly positive or strictly negative norms
# Ouput an orthonormal tetrad with respect to the Minkowski Cartesian metric (+c2,-1,-1,-1)
# Gram-Schmidt orthogonalization, then normalization (c, -1)
# if one or more of the provided vectors is a linear combination of the other vectors,
# the routine will fail with likely outputted NaN
# be aware of the confusion between the names of tetrad and vierbein, due to the definition of a vierbein
# in this routine, two matrices are outputted, tetrad and tetrad_inv, that are in fact the vierbein and inverse vierbein transforms
# according to: chart_vector=tetrad*tetrad_vector and tetrad_vector=tetrad_inv*chart_vector
    (mt_s,mt_si)=Given_Metric(y) # covariant and contravariant metric
    tetrad=fill(BF0,(4,4))
    # Gram-Schmidt construction
    for k=1:4
        for i=1:4
            tetrad[i,k]=SetLinIndVect[i,k]
        end
        for n=2:k
            NUM=BF0
            DEN=BF0
            for p=1:4
                for q=1:4
                    NUM=NUM+mt_s[p,q]*SetLinIndVect[p,k]*tetrad[q,n-1]
                    DEN=DEN+mt_s[p,q]*tetrad[p,n-1]*tetrad[q,n-1]
                end
            end
            for i=1:4
                tetrad[i,k]=tetrad[i,k]-NUM/DEN*tetrad[i,n-1]
            end
        end
    end 
    # orthonormalization
    for k=1:4
        norm2=BF0
        for i=1:4
            for j=1:4
                norm2=norm2+mt_s[i,j]*tetrad[i,k]*tetrad[j,k]
            end
        end
        norm=sqrt(abs(norm2)) # abs(norm2), because normalization means +/- 1, not only 1
        for i=1:4
            tetrad[i,k]=tetrad[i,k]/norm
        end
    end
    for i=1:4
        tetrad[i,1]=tetrad[i,1]*c_light
    end
# verification of validity of tetrad, just a verification of the definition of a tetrad
    tetrad_verif=fill(BF0,(4,4))
    for i=1:4
        for j=1:4
            tetrad_verif[j,i]=BF0
            for p=1:4
                for q=1:4
                    tetrad_verif[j,i]=tetrad_verif[j,i]+mt_s[p,q]*tetrad[p,i]*tetrad[q,j]
                end
            end
        end
    end
    dif2=BF0
    (mt_m,)=Minkowski_Metric_Cartesian(y)
    for i=1:3
        for j=1:4
            dif2=dif2+(tetrad_verif[i,j]-mt_m[i,j])^2
        end
    end
    dif=sqrt(dif2)
    if dif > eps_tetrad # eps_tetrad to be defined by the end-user
        prt3("possible instability in Vierbein routine=",Float64(dif))
    end
    dif=verif_tetrad(tetrad,y,Given_Metric)
    tetradinv=inv(tetrad)
return (tetrad,tetradinv)
end
function verif_tetrad(tetrad::Matrix{BigFloat},y::Vector{BigFloat},Given_Metric::Function)
# verifies if a proposed tetrad is really a tetrad, up to a given accuracy (eps_tetrad)
    (mt_s,)=Given_Metric(y) # covariant and contravariant metric
    (mt_m,)=Minkowski_Metric_Cartesian(y) # covariant and contravariant metric
    # verification of validity of tetrad
    tetrad_v=fill(BF0,(4,4))
    for i=1:4
        for j=1:4
            tetrad_v[j,i]=BF0
            for p=1:4
                for q=1:4
                    tetrad_v[j,i]=tetrad_v[j,i]+mt_s[p,q]*tetrad[p,i]*tetrad[q,j]
                end
            end
        end
    end
    dif2=BF0
    for i=1:3
        for j=1:4
            dif2=dif2+(tetrad_v[i,j]-mt_m[i,j])^2
        end
    end
    dift=sqrt(dif2)
    if dift > eps_tetrad # eps_tetrad to be defined by the end-user
        prt3("possible inaccuracy in tetrad=",Float64(dift))
    end
    tetradinv=inv(tetrad)
    tetrad_v=fill(BF0,(4,4))
    for i=1:4
        for j=1:4
            tetrad_v[j,i]=BF0
            for p=1:4
                for q=1:4
                    tetrad_v[j,i]=tetrad_v[j,i]+mt_m[p,q]*tetradinv[p,i]*tetradinv[q,j]
                end
            end
        end
    end
    dif2=BF0
    for i=1:3
        for j=1:4
            dif2=dif2+(tetrad_v[i,j]-mt_s[i,j])^2
        end
    end
    difti=sqrt(dif2)
    if difti > eps_tetrad # eps_tetrad to be defined by the end-user
        prt3("possible inaccuracy in tetradinv=",Float64(difti))
    end
    return (dift,difti)
end
function FWtetrad_from_y(y::Vector{BigFloat})
    # returns the tetrad/vierbein contained in state vector y, as well as the inverse tetrad
    # according to: chart_vector=tetrad*tetrad_vector and tetrad_vector=tetrad_inv*chart_vector
    tetrad=fill(BF0,(4,4)) 
    for i=1:4 # tetrad vector
        for j=1:4 # vector elements
            tetrad[j,i]=y[4*i+j+6]
        end
    end
    tetradinv=inv(tetrad)
    return (tetrad,tetradinv)
 end
function y_from_FWtetrad(y::Vector{BigFloat},tetrad::Matrix{BigFloat})
# modify the tetrad/vierbein part of y   
# according to: chart_vector=tetrad*tetrad_vector and tetrad_vector=tetrad_inv*chart_vector
    for i=1:4 # tetrad vector
        for j=1:4 # vector elements
            y[4*i+j+6]=tetrad[j,i]
        end
    end
    return y
end
function Cholesky_tetrad_der(y::Vector{BigFloat},Given_Metric::Function)
    # returns: ((tetrad,tetradder=directional derivative of the tetrad along the path,tetradparder=partial derivatives of the tetrad along the path)
    # and same for the inverse of the tetrad
    # Cholesky-like tetrad, with only 10-free parameters, instead of the full 16 free parameters
    # can always be defined for any smooth metric
    tetradinv=fill(BF0,(4,4))
    a=fill(BF0,(4,4))
    M=fill(BF0,(4,4))
    M[1,1]=c_light^2 ; M[2,2]=-BF1 ; M[3,3]=-BF1 ; M[4,4]=-BF1
    (mt_s,mt_si)=Given_Metric(y) # covariant and contravariant metric
    c2=c_light^2
    a[1,1]=sqrt(mt_s[1,1]/c2) ; a[1,2]=mt_s[1,2]/a[1,1]/c2 ; a[1,3]=mt_s[1,3]/a[1,1]/c2 ; a[1,4]=mt_s[1,4]/a[1,1]/c2
    a[2,2]=sqrt(a[1,2]^2*c2-mt_s[2,2]) ; a[2,3]=(a[1,2]*a[1,3]*c2-mt_s[2,3])/a[2,2] ; a[2,4]=(a[1,2]*a[1,4]*c2-mt_s[2,4])/a[2,2]
    a[3,3]=sqrt(a[1,3]^2*c2-a[2,3]^2-mt_s[3,3]) ; a[3,4]=(a[1,3]*a[1,4]*c2-a[2,3]*a[2,4]-mt_s[3,4])/a[3,3]
    # bug corrected
    a[4,4]=sqrt(a[1,4]^2*c2-a[2,4]^2-a[3,4]^2-mt_s[4,4])
    tetradinv[1,1]=a[1,1] ; tetradinv[1,2]=a[1,2] ; tetradinv[1,3]=a[1,3] ; tetradinv[1,4]=a[1,4]
    tetradinv[2,2]=a[2,2] ; tetradinv[2,3]=a[2,3] ; tetradinv[2,4]=a[2,4]
    tetradinv[3,3]=a[3,3] ; tetradinv[3,4]=a[3,4]
    tetradinv[4,4]=a[4,4]
    tetrad=inv(tetradinv)
    b=fill(BF0,(10,10))
    b[1,1] =BF2*c2*a[1,1]
    b[5,2] =BF2*c2*a[1,2]   ; b[5,5] =-BF2*a[2,2]
    b[8,3] =BF2*c2*a[1,3]   ; b[8,6] =-BF2*a[2,3] ; b[8,8] =-BF2*a[3,3]
    b[10,4]=BF2*c2*a[1,4]   ; b[10,7]=-BF2*a[2,4] ; b[10,9]=-BF2*a[3,4] ; b[10,10]=-BF2*a[4,4]
    b[2,1] =c2*a[1,2]       ; b[2,2] =c2*a[1,1]
    b[3,1] =c2*a[1,3]       ; b[3,3] =c2*a[1,1]
    b[4,1] =c2*a[1,4]       ; b[4,4] =c2*a[1,1]
    b[6,2] =c2*a[1,3]       ; b[6,3] =c2*a[1,2]   ; b[6,5]=-a[2,3]      ;  b[6,6]=-a[2,2]
    b[7,2] =c2*a[1,4]       ; b[7,4] =c2*a[1,2]   ; b[7,5]=-a[2,4]      ;  b[7,7]=-a[2,2]
    b[9,3] =c2*a[1,4]       ; b[9,4] =c2*a[1,3]   ; b[9,6]=-a[2,4]      ;  b[9,7]=-a[2,3] ; b[9,8]=-a[3,4] ; b[9,9]=-a[3,3]
    b[10,4]=BF2*c2*a[1,4]   ; b[10,7]=-BF2*a[2,4] ; b[10,9]=-BF2*a[3,4] ;  b[10,10]=-BF2*a[4,4]
    parder_metric_y=Numerical_parder_metric_y(order_parder_metric,y,Given_Metric)
    s_member=fill(BF0,(10,4))
    LINE=0
    for k=1:4
        for l=k:4
            LINE=LINE+1
            for r=1:4
                s_member[LINE,r]=s_member[LINE,r]+parder_metric_y[r,k,l]
            end
        end
    end
    sol=inv(b)*s_member
    tetradparder=fill(BF0,(4,4,4))
    tetradinvparder=fill(BF0,(4,4,4))
    for r=1:4
        tetradinvparder[r,1,1]=sol[1,r] ; tetradinvparder[r,1,2]=sol[2,r] ; tetradinvparder[r,1,3]=sol[3,r] ; tetradinvparder[r,1,4]=sol[4,r]
        tetradinvparder[r,2,2]=sol[5,r] ; tetradinvparder[r,2,3]=sol[6,r] ; tetradinvparder[r,2,4]=sol[7,r]
        tetradinvparder[r,3,3]=sol[8,r] ; tetradinvparder[r,3,4]=sol[9,r]
        tetradinvparder[r,4,4]=sol[10,r]
    end
    tetradinvder=fill(BF0,(4,4))
    for i=1:4
        for j=1:4
            for r=1:4
                tetradinvder[i,j]=tetradinvder[i,j]+tetradinvparder[r,i,j]*y[2*r]
            end
        end
    end
    tetradder=-tetrad*tetradinvder*tetrad
    intermediate=fill(BF0,(4,4))
    for r=1:4
        for i=1:4
            for j=1:4
                intermediate[i,j]=tetradinvparder[r,i,j]
            end
        end
        intermediate=-tetrad*intermediate*tetrad
        for i=1:4
            for j=1:4
                tetradparder[r,i,j]=intermediate[i,j]
            end
        end
    end
    # tetradparder and tetradinvparder can be returned if necessary
    return ((tetrad,tetradder),(tetradinv,tetradinvder))
end
function pieces_of_y(y::Vector{BigFloat})
    # pieces of y (or yp)
    # the velocities are either wrt proper time / integration time / global chart / local tetrad, that does not matter
    position=fill(BF0,4)
    velocity=fill(BF0,4)
    coordinate_time=fill(BF0,2)
    FW_tetrad=fill(BF0,(4,4))
    for i=1:4
        position[i]=y[2*i-1]
        velocity[i]=y[2*i]
    end
    coordinate_time[1]=y[9]
    coordinate_time[2]=y[10]
    for i=1:4 # tetrad vector
        for j=1:4 # elements of each vector tetrad, each column vector = a tetrad vector
            FW_tetrad[j,i]=y[4*i+j+6]
        end
    end
    return (position,velocity,coordinate_time,FW_tetrad)
end
function y_from_pieces(position::Vector{BigFloat},velocity::Vector{BigFloat},coordinate_time::Vector{BigFloat},FW_tetrad::Matrix{BigFloat})
    # y (or yp) from pieces of y (or yp)
    # nsyst should be = length(position)+lenght(velocity)+length(coordinate_time)+length(FW_tetrad)
    # the velocities are either wrt proper time / integration time / global chart / local tetrad, that does not matter
    y=fill(BF0,nsyst)
    for i=1:4
        y[2*i-1]=position[i]
        y[2*i]=velocity[i]
    end
    y[9]=coordinate_time[1]
    y[10]=coordinate_time[2]
    for i=1:4 # tetrad vector
        for j=1:4 # elements of each vector tetrad, each column vector = a tetrad vector
            y[4*i+j+6]=FW_tetrad[j,i]
        end
    end
    return y
end
function velocity_from_y(y::Vector{BigFloat})
    # extract the velocity vector from y
    velocity=fill(BF0,4)
    for p=1:4
        velocity[p]=y[2*p]
    end
    return velocity
    end
function chart_vector_from_tetrad_vector(local_vector::Vector{BigFloat},y::Vector{BigFloat})
    # from vector components wrt to tetrad frame to vector components wrt the local holonomic frame
    (tetrad,tetradinv)=FWtetrad_from_y(y)
    chart_vector=fill(BF0,4)
    for p=1:4
        for q=1:4
            chart_vector[p]=chart_vector[p]+tetrad[p,q]*local_vector[q]
        end
    end
    return chart_vector
end
function tetrad_vector_from_chart_vector(chart_vector::Vector{BigFloat},y::Vector{BigFloat})
    # from vector components wrt the local holonomic frame to vector components wrt to tetrad frame
    (tetrad,tetradinv)=FWtetrad_from_y(y)
    tetrad_vector=fill(BF0,4)
    for p=1:4
        for q=1:4
            tetrad_vector[p]=tetrad_vector[p]+tetradinv[p,q]*chart_vector[q]
        end
    end
    return tetrad_vector
end
function ychol_from_ychart(y::Vector{BigFloat},Given_Metric)
    # position given in chart coordinates,
    # velocity given in chart coordinates (chart tetrad to be totally correct), any time of integration
    (position_in_chart,velocity_in_chart,dummy1,dummy2)=pieces_of_y(y)
    ((tetradchol,tetradcholder),(tetradcholinv,tetradcholinvder))=Cholesky_tetrad_der(y,Given_Metric)
    velocity_in_chol_tetrad=tetradcholinv*velocity_in_chart # from velocity in chart to local Cholesky tetrad
    ychol=y_from_pieces(position_in_chart,velocity_in_chol_tetrad,dummy1,dummy2)
    return ychol
end
function ychart_from_ychol(ychol::Vector{BigFloat},Given_Metric)
    # position given in chart coordinates,
    # velocity given in Cholesky tetrad, any choice of time of integration
    (position_in_chart,velocity_in_chol_tetrad,dummy1,dummy2)=pieces_of_y(ychol)
    # only position (so in chart coordinates) is used inside Cholesky_tetrad_der, so we can input ychol in Cholesky_tetrad_der
    ((tetradchol,tetradcholder),(tetradcholinv,tetradcholinvder))=Cholesky_tetrad_der(ychol,Given_Metric)
    velocity_in_chart=tetradchol*velocity_in_chol_tetrad
    chart_y=y_from_pieces(position_in_chart,velocity_in_chart,dummy1,dummy2)
    return chart_y
end
function ypchol_from_ypchart(y::Vector{BigFloat},yp::Vector{BigFloat},Given_Metric)
    # velocity in chart coordinates
    # acceleration in chart coordinates
    # (position_in_chart,velocity_in_chart0,dmmy1,dumm2)=pieces_of_y(y::Vector{BigFloat})
    (velocity_in_chart,acc_in_chart,dummy1,dummy2)=pieces_of_y(yp::Vector{BigFloat})
    # caution: tetradcholinvder is wrt the time used for the time derivative
    ((dummy3,dummy4),(tetradcholinv,tetradcholinvder))=Cholesky_tetrad_der(y,Given_Metric)
    acc_in_chol_tetrad=tetradcholinvder*velocity_in_chart+tetradcholinv*acc_in_chart
    ypchol=y_from_pieces(velocity_in_chart,acc_in_chol_tetrad,dummy1,dummy2)
    return ypchol
end
function ypchart_from_ypchol(ychol::Vector{BigFloat},ypchol::Vector{BigFloat},Given_Metric)
    # (position_in_chart,velocity_in_chart0,dmmy1,dumm2)=pieces_of_y(ychol::Vector{BigFloat})
    yt=ychart_from_ychol(ychol::Vector{BigFloat},Given_Metric)
    (velocity_in_chart,acc_in_chol_tetrad,dummy1,dummy2)=pieces_of_y(ypchol::Vector{BigFloat})
    ((tetradchol,tetradcholder),(tetradcholinv,dummy3))=Cholesky_tetrad_der(yt,Given_Metric)
    # caution: tetradcholder is wrt the time used for the time derivative
    acc_in_chart=tetradcholder*tetradcholinv*velocity_in_chart+tetradchol*acc_in_chol_tetrad
    ypchart=y_from_pieces(velocity_in_chart,acc_in_chart,dummy1,dummy2)
    return ypchart
end
function Non_Grav_Grad(y::Vector{BigFloat})
    # non-gravitational gradient (a tensor) of the energy-momentum tensor, as only a function of 4-chart coordinates
    # only the part that is orthogonal to the worldline is "seen" by the spacecraft
    NGG=fill(BF0,4) # contravariant tensor wrt proper time
    # here an example for Parker Solar probe with actual values (radiation pressure from the Sun), see SoftwareX paper
    NGG[1]=BF0
    # const mass_PSP=parse(BigFloat,"655.0") # mass in kg, PSP
    # const shield_area=parse(BigFloat,"4.0") # Sun shield area in m^2, PSP
    # const shiel_reflection_coef=parse(BigFloat,"1.8") # shield reflection coef, shield always facing the Sun at periapsis
    # radiation along the r polar chart coordinate here
    NGG[2]=solar_radiation_pressure*shield_area*shield_reflection_coef*(AU/y[3])^2/mass_PSP
    NGG[3]=BF0
    NGG[4]=BF0
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # NGG=fill(BF0,4) # put zero for the case of non-gravitational forces, so no NG forces if this line is active,
    # used for modeling the Lense-Thirring effect in the SoftwareX paper
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return NGG
end
function Non_Grav_Bus(y::Vector{BigFloat})
    # non-gravitational forces measured on the spacecraft bus 
    # the first component MUST BE, from physics, ALWAYS zero (see SoftwareX paper)
    NGB=fill(BF0,4)
    NGB[1]=BF0 # DO NOT CHANGE IT, THIS COMPONENT IS ALWAYS ZERO BY ***DEFINITION*** IN THE S/C BUS
    NGB[2]=BF0 # as you want (model or measurement by S/C onboard accelerometers)
    NGB[3]=BF0 # as you want (idem)
    NGB[4]=BF0 # as you want (idem)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NGB=fill(BF0,4) # put to zero for the test case of SotwareX, so no NG forces here if this line is active
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return NGB
end
function Eqs_Motion_chol(evpar::BigFloat,ychol::Vector{BigFloat},Given_Metric::Function,Non_Grav_Grad::Function,Non_Grav_Bus::Function)
    # the spherical coordinates integrated here are non-observable quantities
    # K contravariant gradient of energy-momentum tensor
    # reminder: the spherical coordinates integrated here are non-observable quantities
    # contravariant coordinates and velocities in spatial polar coordinates here
    # t=y[1]=coordinate time (clock at infinity) ; r=y[3]=pseudo-radius
    # phi=y[5]=pseudo-longitude ; theta=y[7]=pseudo-colatitude,
    # tau is proper time of S/C
    nsyst=length(ychol)
    ypchol=fill(BF0,nsyst)
    yp=fill(BF0,nsyst)
    # going back to y from ychol for the computations in Eqs_Motion_Chol !!!!!!!!!!!!!!!!!!!
    ychart=ychart_from_ychol(ychol,Given_Metric) # so in the integration time
    chart_four_velocity=velocity_from_y(ychart) # here, it is 4-velocity,i.e. wrt to proper time
    chol_velocity=velocity_from_y(ychol) # here it is velocity from ychol, wrt to the time of integration
    ((tetradchol,tetradcholder),(tetradcholinv,tetradcholinvder))=Cholesky_tetrad_der(ychart,Given_Metric) # so here wrt integration time
    dtaudevpar_external=evpar2tau_external(time_option,evpar,ychart,Given_Metric) # this give from an external routine a priori the values of tau and its derivative wrt evpar
    ##-------------------------------------------------------------------------------------------------
    dtaudevpar_eqm=dtaudevpar_external
    ##-------------------------------------------------------------------------------------------------
    (mt_s,mt_si)=Given_Metric(ychart) # covariant and contravariant metric
    NGG=Non_Grav_Grad(ychart) # covariant gradient (provided) of the energy-momentum tensor
    NGB=Non_Grav_Bus(ychart) # Non-Gravitational forces measured wrt the spacecraft
    # order_parder_metric is defined in main (see user's guide)
    BigGamma=Numerical_Christoffel_Symbols(order_parder_metric,ychart,Given_Metric)
    F=fill(BF0,4) # contravariant 4-force, orthogonal to the wordline
    for i=1:4
        for j=1:4
            F[i]=F[i]+NGG[j]*(mt_si[j,i]-chart_four_velocity[j]*chart_four_velocity[i]/c_light^2)
        end
    end
    ypBG=fill(BF0,4)
    yp_tot=fill(BF0,4)
    yp_pt=fill(BF0,nsyst)
    for i=1:4
        for j=1:4
            for k=1:4
                ypBG[i]=ypBG[i]+BigGamma[i,j,k]*chart_four_velocity[j]*chart_four_velocity[k] # Big Gamma gravity acceleration
            end
        end
    end
    for i=1:4
        yp_tot[i]=ypBG[i]-F[i] # sum of Big Gamma gravity acceleration + non-gravitational forces (with a minus sign)
        yp_pt[2*i]=yp_tot[i]
    end
    ypchol_vel=-dtaudevpar_eqm*tetradcholinv*(tetradcholder*chol_velocity+yp_tot) # as ypBG and F are evaluated wrt proper time, tetradcholder wet to integration time (ychart)
    for i=1:4
        ypchol[2*i-1]=dtaudevpar_eqm*ychart[2*i] # careful, this is because we need to have the chart coordinates with respect to integration time
        ypchol[2*i]=ypchol_vel[i]
    end
    ypchol[9]=dtaudevpar_eqm # normally, dtaudevpar_external, so will give proper time by integration
    ypchol[10]=BF1 # integration time (elapsed)
# ------------------------------------------------computation of the transported tetrad---------------------------------------------------------------------------
# parallel transport wrt integration time of the spatial part of the tetrad
# caution: this will only work as a Fermi-walker transport if and only if the first vector of the tetrad is the 4-velocity
# otherwise the transported tetrad will rotate / wrt distant stars
    for p=1:4
        for i=1:4
            yp[4*p+i+6]=BF0
            for j=1:4
                for k=1:4
                    yp[4*p+i+6]=yp[4*p+i+6]-BigGamma[i,j,k]*ychart[4*p+j+6]*ychart[2*k]
                end
            end
        end
    end
# Khl = correction to the parallel transport to get the Fermi-Walker transport
    Khh=fill(BF0,(4,4))
    for p=1:4
        for q=1:4
            Khh[p,q]=(F[p]*ychart[2*q]-F[q]*ychart[2*p])/c_light^2
        end
    end
    Khl=fill(BF0,(4,4))
    for p=1:4
        for q=1:4
            for l=1:4
                Khl[p,q]=Khl[p,q]+mt_s[l,q]*Khh[p,l]
            end
        end
    end
# application of the Khl correction to get the Fermi-Walker transport
    for p=1:4
        for i=1:4
            for l=1:4
                yp[4*p+i+6]=yp[4*p+i+6]+Khl[i,l]*ychart[4*p+l+6]
            end
        end
    end
    for i=11:26
        ypchol[i]=dtaudevpar_eqm*yp[i]
    end
    return ypchol
end
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
function rki(t0::BigFloat,t1::BigFloat,y0::Vector{BigFloat},Given_Metric::Function,Eqs_Motion::Function,Non_Grav_Grad::Function,Non_Grav_Bus::Function)
    # apply to Eqs_Motion or Eqs_Motion_chol indifferently
    # Gauss 5-stages implicit symplectic integrator
    # see Butcher's paper (page 57, Math. Comput. 18(85), 50-64, 1964, ISSN 0025-5718) 
        (stages,)=size(a_B) ; nsyst=length(y0)
        # the k_B's are initialized to zero, pretty rough, but it works,
        # the radius of convergence seems to be large, to be asserted (future work)
        k_B=fill(BF0,(stages,nsyst))
        h_B=t1-t0
        # the implicit method is implemented by a pretty rough fixed point method,
        # with an a priori fixed number of loops (itmax_rki) 
        # itmax_rki=10 for Parker Solar Probe, defined in main
        # to be implemented: better initialization, convergence criteria following
        # the Hairer et al. 2002 book, p 275, future work (see SoftwareX paper)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
        for it=1:itmax_rki
            for i=1:stages
                s_B=fill(BF0,(nsyst,1))
                for j=1:stages
                    for n=1:nsyst
                        s_B[n]=s_B[n]+a_B[i,j]*k_B[j,n]
                    end
                end
                yf=Eqs_Motion(t0+h_B*c_B[i,1],vec(y0+h_B*s_B),Given_Metric,Non_Grav_Grad,Non_Grav_Bus)
                for n=1:nsyst
                    k_B[i,n]=yf[n]
                end
            end
        end 
        dy=(b_B*k_B)'
        y1=vec(y0+h_B*dy)
        return y1
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
function main(evpar0::BigFloat,ychol0::Vector{BigFloat},Given_Metric::Function)
    nsyst=length(ychol0)
    y0=ychart_from_ychol(ychol0,Given_Metric)
    dtaudevpar_internal=BF0 # do not remove, needed for scope of variable 
    y1=fill(BF0,nsyst) # do not remove, otherwise this parameter is only defined inside the next loop !!!!
    dtaudevpar=BF0 # do not remove, otherwise this parameter is only defined inside the next loop !!!!
    ev_par0=evpar0 # just for integration # start of time step
    prt3("entering main loop at=",now())
        open("SCephemeris_$(timetag).txt","w") do f
            for i=1:N_integration_steps
                for j=1:nsyst
                    write(f,"$(j) "*"$(y0[j]) "*"\n")
                end     
                ev_par1=ev_par0+evpar_step # end of time step
                ychol1=rki(ev_par0,ev_par1,ychol0,Given_Metric,Eqs_Motion_chol,Non_Grav_Grad,Non_Grav_Bus) # integration in integration time, quantities given at time=evpar0
                y1=ychart_from_ychol(ychol1,Given_Metric)
                dtaudevpar_internal=der1_tau_internal(y1,Given_Metric)
                ev_par0=ev_par1 ; y0=deepcopy(y1) ; ychol0=deepcopy(ychol1)    
            end # i loop  
            for j=1:nsyst
                write(f,"$(j) "*"$(y0[j]) "*"\n")
            end
        end # open
        prt3("exiting main loop at=",now())
     return ev_par0,ychol0,dtaudevpar_internal # tuple of values returned by main()
end
# --------------------------------------------------------------root of module-------------------------------------------------------------------------------
#-------------------------------------the constants defined here propagate throught the entire code----------------------------------------------------------
#---------------------------------------------------Windows 11 implementation--------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
const BF0=BigFloat(0)       ; const BF1=BigFloat(1)       ; const BF6=BigFloat(6)
const BF2=BigFloat(2)       ; const BF8=BigFloat(8)       ; const BF11=BigFloat(11)
const BF13=BigFloat(13)     ; const BF23=BigFloat(23)
const BF32=BigFloat(32)     ; const BF35=BigFloat(35)     ; const BF59=BigFloat(59)
const BF63=BigFloat(63)     ; const BF64=BigFloat(64)     ; const BF70=BigFloat(70)
const BF225=BigFloat(225)   ; const BF308=BigFloat(308)   ; const BF322=BigFloat(322)
const BF405=BigFloat(405)   ; const BF452=BigFloat(452)   ; const BF960=BigFloat(960)
const BF1080=BigFloat(1080) ; const BF3240=BigFloat(3240) ; const BF3600=BigFloat(3600)
const BF05=BF1/BF2
const BF3=BigFloat(3) ; const BF4=BigFloat(4) ; const BF5=BigFloat(5) ; const BF12=BigFloat(12) ; const BF20=BigFloat(20)
const BF60=BigFloat(60) ; const BF105=BigFloat(105) ; const BF280=BigFloat(280)	
const BF180=BigFloat(180) ; const BF1000=BigFloat(1000)
# Gauss 5-stages implicit symplectic integrator
# see Butcher's paper (page 57, Math. Comput. 18(85), 50-64, 1964, ISSN 0025-5718) 
const om1=(BF322-BF13*sqrt(BF70))/BF3600        ; const op1=(BF322+BF13*sqrt(BF70))/BF3600
const om2=BF05*sqrt((BF35+BF2*sqrt(BF70))/BF63) ; const op2=BF05*sqrt((BF35-BF2*sqrt(BF70))/BF63)
const om3=om2*(BF452+BF59*sqrt(BF70))/BF3240    ; const op3=op2*(BF452-BF59*sqrt(BF70))/BF3240
const om4=om2*(BF64+BF11*sqrt(BF70))/BF1080     ; const op4=op2*(BF64-BF11*sqrt(BF70))/BF1080
const om5=BF8*om2*(BF23-sqrt(BF70))/BF405       ; const op5=BF8*op2*(BF23+sqrt(BF70))/BF405
const om6=om2-BF2*om3-om5                       ; const op6=op2-BF2*op3-op5
const om7=om2*(BF308-BF23*sqrt(BF70))/BF960     ; const op7=op2*(BF308+BF23*sqrt(BF70))/BF960
const a_B=[om1          op1-om3+op4  BF32/BF225-om5  op1-om3-op4  om1-om6      ; # 2D matrix
           om1-op3+om4  op1          BF32/BF225-op5  op1-op6      om1-op3-om4  ;
           om1+om7      op1+op7      BF32/BF225      op1-op7      om1-om7      ;
           om1+op3+om4  op1+op6      BF32/BF225+op5  op1          om1+op3-om4  ;
           om1+om6      op1+om3+op4  BF32/BF225+om5  op1+om3-op4  om1          ]
const b_B=[BF2*om1  BF2*op1  BF64/BF225  BF2*op1  BF2*om1] # line-vector (1-line matrix, adressing[1,*] or [*]), size(b)=(1,*)
const c_B=[BF05-om2 , BF05-op2 , BF05 , BF05+op2 , BF05+om2] # column-vector 
# global scope of module
using Dates # just to get access to the computer clock
using PythonPlot # this is the plotting library, can be replaced with minimal work
# using SPICE # to be activated to use de430, with a dedicated package
td="$(now())"
timetag=td[1:4]*td[6:7]*td[9:10]*td[12:13]*td[15:16]*td[18:19] # to tag the figures with time 
mkdir("EXGR_"*timetag)
cd("EXGR_"*timetag)
runlogfile=open("runlog_$(timetag).txt","w")
prt3("(G)eneral (R)elativity (A)ccelerometer-based (P)ropagation (E)nvironment, version April 12, 2026")
prt3("timetag=",timetag,"==> all the files related to this run are labeled with this timetag") # important, gives the link between this run and the labels of the figures and logfile
prt3("and outputed in the working directory thereafter")
working_directory=pwd() # to know in which directory the computation is taking place
prt3("working directory=",working_directory)
dumpcode(timetag) # copy the running code into a file named codefile_timetag.txt
t_local=now()
prt3("computer clock at start of run=",t_local)
#--------------------------------------------------------------------------------------------------------------------------------------------------------
const nsyst=26 # number of systems to integrate (8: eqs of motion, 2: proper time, 16: tetrad)
#--------------------------------------------------------------------------------------------------------------------------------------------------------
setprecision(BigFloat,237) # 53 for 64 bits IEEE754, 113 for 128 bits, 237 for 256 bits IEEE754, using the GNU-MPFR library
precision_bits=precision(BigFloat)
prt3("BigFloats bits for MPFR library=",precision_bits)
prt3("control of accuracy of MPFR library with respect to external value of PI ==>")
pi_run=string(one(BigFloat)*π)
prt3(pi_run,length(pi_run)-1,"digits for this run from MPFR library")
pi_ref="3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342"
prt3(pi_ref,length(pi_ref)-1,"digits of pi for control") # to check accuracy of the computation of BigFloats
# NEVER do this type of stuff in the code: BigFloat(2.1) or 2.1*BF1, except for pi that is autoadaptative, for Julia versions < 1.9
# you will get this ===> 2.100000000000000088817841970012523233890533447265625000000000000000000000000000, please use Julia >= 1.9
#--------------------------------------------------------------------------------------------------------------------------------------------------
# the implicit Runge-Kutta symplectic method used here is implemented by a pretty rough fixed point method,
# with an a priori fixed number of loops (itmax_rki)
# The Runge-Kutta integrator used here is totally symplectic, in the sense that it conserves the norm of the 4-velocity
const itmax_rki=10 # should be usually larger than 6, 10 is needed for Parker Solar Probe in a highly elliptical orbit
prt3("number of iterations for implicit integrator=",itmax_rki)
# the integrator can be changed at will by the end user, but this one is symplectic
order_parder_metric=4 # this defines the precision for the computation (2, 4, 6 or 8) of the numerical derivation of the metric / Christoffel (connection) symbols throughout the code
prt3("central difference order for Christoffel symbols=",order_parder_metric)
# epsmin_CSN, epsabs are used in Numerical_parder_metric_y, these values must be checked to be sure that they are OK for the case being studied
const epsabs_CSN=parse(BigFloat,"1.e-12")  # cte for numerical derivation of Christoffel (connection) symbols
const epsmin_CSN=parse(BigFloat,"1.e-50")  # cte for numerical derivation of Christoffel (connection) symbols
prt3("first constant for numerical derivation of Christoffel symbols=",Float64(epsabs_CSN))
prt3("second constant for numerical derivation of Christoffel symbols=",Float64(epsmin_CSN))
const eps_tetrad=parse(BigFloat,"1.e-9") # to check validity of the outputted tetrad in verif_tetrad, this will not stop the computation
prt3("threshold for the numerical test of the validity of Fermi-Walker tetrad=",Float64(eps_tetrad))
# this is the choice of the integration time as: (1: proper time, 2: function of proper time, 3: coordinate time, 4: function of coordinate time)
const time_option=3 # the integration time can be either proper time or coordinate time, or a function of them (see routine evpar2tau_external)
prt3("integration time option=",time_option) # =1, proper time, =2, user defined time wrt proper time, =3, coordinate time
prt3("1 ==> proper time, 2 ==> user defined time wrt proper time, 3 ==> coordinate time, 4 ==> user defined time wrt coordinate time")
const evpar_step=parse(BigFloat,"120.00") # integration time step in seconds for the integration process
prt3("integration time step=",Float64(evpar_step),"seconds")
const evpar0=BF0 # integration time start value for the integration process
prt3("integration time at start of integration=",Float64(evpar0),"seconds")
N_integration_steps=250_000 # number of time steps for the integration process
prt3("number of loops for integration process=",N_integration_steps)
prt3("duration period of probe orbit in this run=",Float64(N_integration_steps*evpar_step),"seconds of integration time")
c_light_nature=parse(BigFloat,"299792.458") # light velocity in km, fixed by IAU definition
#prt3("c_light_in_nature=",Float64(c_light_nature),"km s-1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
const c_light=c_light_nature # natural constant value, was divided by 10 for tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prt3("c light=",Float64(c_light),"km s^-1")
const eps_c2=BF1/c_light^2 # used instead of 1/c2, eps_c2=0 Minkowskian metric
const Kerr_Metric_Polar_name=Give_Function_Name(Kerr_Metric_Polar) # do not remove this line
const G=parse(BigFloat,"6.6743e-20") # cte of gravitation (unit: km3/kg/s2)
#prt3("constant of gravitation G=",Float64(G),"km^3 kg^-1 s^-2")
AU=parse(BigFloat,"149597870.700") # Astronomical Unit in km
#prt3("one astronomical unit=",Float64(AU),"km")
#Given_Metric=Give_Function_Name(Schwarzschild_Metric_Polar) # another case of metric for tests
const mu_s=parse(BigFloat,"0.132712440018e12") # gm of Sun in km3/s2
#prt3("GM of Sun=",mu_s,"km^3 s^-2")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
const mu=mu_s # gm of star 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#mu=BF0 # to go back to Minkowski space-time in spherical coordinates, if needed
prt3("GM of star=",Float64(mu),"km^3 s^-2")
const Mstar=mu/G
prt3("mass of star=",Float64(Mstar),"kg")
const rs=eps_c2*BF2*mu # Schwarzschild radius km
prt3("Schwarzschild radius of star=",Float64(rs),"km") 
# angular momentum of the Sun from helioseismology: 2.02 +/- 0.04 10^48 g cm^2 s^-1, Mauro et al. (2000), ADS 2000ASPC..198..353D
const Given_Metric=Give_Function_Name(Kerr_Metric_Polar)
prt3("Given_Metric=",Given_Metric)
if Given_Metric == Kerr_Metric_Polar_name
   prt3("initialization of Kerr metric")
   const Jsun=parse(BigFloat,"2.02e35") # Angular momentum of the Sun kg km^2 s^-1 # value=1.92 1o^41 kg m^2 s^-1
   prt3("angular momentum of the Sun=",Float64(Jsun),"km^2 kg s^-1")
   const Jstar=Jsun
   prt3("angular momentum of star=",Float64(Jstar),"km^2 kg s^-1")
   const as=Jstar/Mstar/c_light # this is the constant in length units entering the Kerr's metric
   # const as=BF0  # if needed, to go back to a Schwarschild metric
   prt3("scaled angular momentum of star=",Float64(as),"km") # as (in km) is used in the definition of the Kerr metric, it is the angular momentum per unit of mass scaled wrt speed of light
   const event_horizon=rs/BF2+sqrt(rs^2/BF4-as^2) # even horizon of the ergosphere
   prt3("event horizon of star=",Float64(event_horizon),"km")
else
   const as=BF0
end
const solar_constant=parse(BigFloat,"1361.0") # solar constant at 1 AU in Watt/m^2, exact value between 1360. and 1367. depending on references, fluctuates with the solar cycle of 11 years
solar_radiation_pressure=solar_constant/c_light/BF1000 # solar pressure in Newtons/m^2
prt3("solar_radiation_pressure=",Float64(solar_radiation_pressure*BF1000*BF1000),"micro-N/m^2")
# initial and running state vectors: eqs of motion 8; proper time 2; tetrad 16
y0_pt=fill(BF0,nsyst) # contravariant vector, initial state vector must be *wrt proper time* in chart coordinates
# all the stuff between the ~~~~~~~~~~~~~~~~~~~~~~ lines is the initialisation vector for the Parker Solar probe, can be changed at will 
# ------------------------------------------------------------------------------------------------------------------------------------------
#=
Solar shield reflection coefficient 1.8
Solar panels reflection coefficient 1.38
Solar shield’s area 4 m^2
Solar panel’s area 1.6 m^2
Parker Solar Probe’s mass 655 kg
Intensity of Solar radiation at 1 AU: 1367 W =m2
=#
const mass_PSP=parse(BigFloat,"655.0") # mass in kg, PSP
const shield_area=parse(BigFloat,"4.0") # Sun shield area in m^2, PSP
const shield_reflection_coef=parse(BigFloat,"1.8") # shield reflection coef, shield always facing the Sun at periapsis
#=
y0_pt is the initial state vector per se, given ****wrt to proper time**** of the spacecraft
these values are given with respect to the Kerr polar metric provided in this code
y0[1]=t0 # coordinate time initial value # must be coordinate time, ****cannot be changed****
y0[2]=dtdtau0 # initial rate of coordinate time wrt proper time
# for the space part, you can use any sysyem of coordinates you want, provided that you are consistent throughout the software
y0[3]=r0 # radius distance initial value
y0[4]=drdtau0 # initial r rate wrt proper time
y0[5]=phi0 # longitude initial value
y0[6]=dphidtau0 # initial phi rate wrt proper time
y0[7]=theta0 # colatitude initial value
y0[8]=dthetadtau0 # initial theta rate proper time
y0[9]=tau0 # proper time initial value
y0[10]=evpar0 # integration time initial value
=#
# BODY_GM=init_de430() # to be activated in order to use de430 ephemeris
#----------------------------------initial vector for PSP UT Time: 2025-06-09T00:00:00.00 --------------------------------------------
# pure heliocenric orbit
y0_pt[1]=parse(BigFloat,"8.026992691847092000000000000000000000000000000000000000000000000000000013e+08") # J2000 initial epoch
y0_pt[2]=parse(BigFloat,"1.000000036065583982352343088267053815632478870312984159711699273598484404")
y0_pt[3]=parse(BigFloat,"6.05555029296586368735956299733688278740147641652955018507419361797892692e+07")
y0_pt[4]=parse(BigFloat,"-39.20853470649714613562963656044536896149439334319005675429975900686120154")
y0_pt[5]=parse(BigFloat,"-0.1051484868287112526396508187748939473128207329184157695124121577662568626")
y0_pt[6]=parse(BigFloat,"3.269982323577800306566017661024433703896072724496301294340993590072270686e-07")
y0_pt[7]=parse(BigFloat,"1.570796111324400805115005611947296536407776016865991878896744220327295114")
y0_pt[8]=parse(BigFloat,"2.154706430982307446885933453094922768491406204024097910843647327955235975e-07")
#--------------------------------------------------------------------------------------------------------------------------------------
y0_pt_Newton=fill(BF0,nsyst) # velocities wrt to coordinate time
for i=1:4
    y0_pt_Newton[2*i]=y0_pt[2*i]/y0_pt[2]
end
# dtaudevpar_internal must be equal to unity wrt to numerical precision for y0_pt[2,4,6,8] to be acceptable as an initial vector
# if it is not, this means that the initial values refer to another metric than the metric used for the runs
dtaudevpar_internal=der1_tau_internal(y0_pt,Kerr_Metric_Polar) # must be equal to unity for initial values wrt proper time to be acceptable as an initial vector
prt3("the value at the following line must be equal 1 wrt numerical precision for the initial velocities of the s/c to be physically acceptable:")
prt3("partial derivative of proper time wrt proper time at start (initial vector)=",dtaudevpar_internal)
# the routine evpar2tau_external, given by the end user, links the proper time and its derivatives wrt the integration time
const tau0=BF0 # proper time of the integration arc, starts at zero here, passed as a global variable to routine evpar2tau_external
dtaudevpar_external=evpar2tau_external(time_option,evpar0,y0_pt,Given_Metric)
prt3("partial derivative of proper time wrt integration time at start (imposed)=",dtaudevpar_external)
##==========================================================================================================================
y0_pt[9]=tau0 # proper time
y0_pt[10]=evpar0 # integration time
##==========================================================================================================================
# setup of a comoving tetrad (Fermi-Walker (FW) transport of the orientation of the spacecraft)------------------------------------------------------------------------------------------
# list of independant 4-vectors, the first vector must be proportional to the 4-velocity to have a Fermi-Walker transport, for inputting the Gram-Schmidt orthogonalization process
SetLinIndVect=fill(BF0,(4,4)) 
SetLinIndVect[1,1]=y0_pt[2] ; SetLinIndVect[2,1]=y0_pt[4] ; SetLinIndVect[3,1]=y0_pt[6] ; SetLinIndVect[4,1]=y0_pt[8]
SetLinIndVect[1,2]=BF0      ; SetLinIndVect[2,2]=BF1      ; SetLinIndVect[3,2]=BF0      ; SetLinIndVect[4,2]=BF0
SetLinIndVect[1,3]=BF0      ; SetLinIndVect[2,3]=BF0      ; SetLinIndVect[3,3]=BF1      ; SetLinIndVect[4,3]=BF0
SetLinIndVect[1,4]=BF0      ; SetLinIndVect[2,4]=BF0      ; SetLinIndVect[3,4]=BF0      ; SetLinIndVect[4,4]=BF1
(tetrad,tetradinv)=Vierbein(SetLinIndVect,Given_Metric,y0_pt) # ==> output the local co-moving FW tetrad (at initial coordinates) and its inverse
y0_pt=y_from_FWtetrad(y0_pt,tetrad) # input tetrad in y0_pt
(tetrad,)=FWtetrad_from_y(y0_pt) # just to verify
(dift,difti)=verif_tetrad(tetrad,y0_pt,Given_Metric)
prt3("validation of Fermi-Walker (FW) tetrad at start=",Float64(log10(dift)),Float64(log10(difti)),"(log10 of precision)") # verifies if the outputted candidate tetrad satisfies the definition of a tetrad
chart_4_velocity=velocity_from_y(y0_pt)
(mt_s,mt_si)=Given_Metric(y0_pt)
invariant0=mt_s*chart_4_velocity
ychol0_pt=ychol_from_ychart(y0_pt,Given_Metric) # integration in Cholesky tetrad (not the Fermi-Walker tetrad)
transmit_frequency_Doppler=parse(BigFloat,"8.0e9") # typical frequency for Doppler computation (normally X-band) in Hertz
prt3("transmit frequency for Doppler=",Float64(transmit_frequency_Doppler),"Hertz")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
t_start=now()
prt3("integration start at=",t_start)
((tetradchol,tetradcholder),(tetradcholinv,tetradcholinvder))=Cholesky_tetrad_der(y0_pt,Given_Metric) # this is the tetrad for implementing symplecticity
(dift,difti)=verif_tetrad(tetradchol,y0_pt,Given_Metric)
prt3("validation Cholesky tetrad at start=",Float64(log10(dift)),Float64(log10(difti)),"(log10 of precision)") # tetrad and inverse tetrad precisions
evpar_start=evpar0
evpar_end,ychol_end,dtaudevpar_internal_end=main(evpar0,ychol0_pt,Given_Metric) # ----------------------------------MAIN RUN-------------------------------------
y_end=ychart_from_ychol(ychol_end,Given_Metric)
(tetrad_end,)=FWtetrad_from_y(y_end)
(dift,difti)=verif_tetrad(tetrad_end,y_end,Given_Metric)
dift=log10(dift) ; difti=log10(difti)
prt3("validation Fermi-Walker tetrad at end=",Float64(dift),Float64(difti),"(log10 of precision)")
chart_4_velocity=velocity_from_y(y_end)
(mt_s,mt_si)=Given_Metric(y_end) # covariant and contravariant metric
invariant=mt_s*chart_4_velocity
prt3("covariant momenta at start and end of integration ==> next 5 lines")
prt3("first and third components are constants of motion for a pure Kerr spacetime with no non-gravitational forces")
for i=1:4
    prt3(i,invariant[i],invariant0[i]) # first and third components are constants of motion for a pure Kerr spacetime with no NG forces
end
local_velocity=tetrad_vector_from_chart_vector(chart_4_velocity,y_end) # 4-velocity in tetrad frame
norm2=Minkowski_norm2(local_velocity)
prt3("control 4-velocity norm in Fermi-Walker tetrad, should be close to 1 (wrt c_light_software^2)=",norm2/c_light^2)
t_end=now()
prt3("computation end=",t_end)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
prt3("coordinate time at start and end=",Float64(y0_pt[1]),"seconds,",Float64(y_end[1]),"seconds,","dif=",Float64(y_end[1]-y0_pt[1]))
prt3("proper time at start and end=",Float64(y0_pt[9]),"seconds,",Float64(y_end[9]),"seconds,","dif=",Float64(y_end[9]-y0_pt[9]))
prt3("integration time at start and end=",Float64(y0_pt[10]),"seconds,",Float64(y_end[10]),"seconds,","dif=",Float64(y_end[10]-y0_pt[10]))
for i=1:nsyst
    prt3("initial/final state vector (chart)=",i,y0_pt[i],y_end[i])
end
dtaudevpar_internal_end=der1_tau_internal(y_end,Given_Metric) # should be close to 1.
dtaudevpar_external_end=evpar2tau_external(time_option,evpar_end,y_end,Given_Metric)
prt3("partial derivative of proper time wrt proper time at end (propagated)=",dtaudevpar_internal_end)
#write(runlogfile,prtstr("partial derivative of proper time wrt proper time at end (propagated)=",dtaudevpar_internal_end))
prt3("partial derivative of proper time wrt integration time at end (imposed)=",dtaudevpar_external_end)
#write(runlogfile,prtstr("partial derivative of proper time wrt integration time at end (imposed)=",dtaudevpar_external_end))
prt3("elapsed time for integration=",t_end-t_start)
#------------------------------------------------plots section-----------------------------------------------------------------------------
time_coord_plot=Float64[]
r_plot=Float64[]
dr_plot=Float64[]
phi_plot=Float64[]
theta_plot=Float64[]
tau_plot=Float64[]
time_coord_plot=Float64[]
dif_tau_timecoord_plot=Float64[]
Doppler_plot=Float64[]
dif_tetrad_plot=Float64[]
X_Plot=Float64[]
Y_Plot=Float64[]
Z_Plot=Float64[]
four_velocity_check_plot=Float64[]
y1=fill(BF0,nsyst)
evpar_count=BF0
open("SCephemeris_$(timetag).txt","r") do f # reading the spacecraft ephemeris for plotting
    while !eof(f) # operations for the plot in Float64
        dif_tetrad=BF0
        for j=1:nsyst
            line=readline(f)
            line=rstrip(lstrip(line)) # remove heading and trailing blanks
            fields=split(line) # split(line,r" {1,}") regex expression
            y1[j]=parse(BigFloat,fields[2])
        end
        (tetrad,)=FWtetrad_from_y(y1)
        (dif_tetrad,)=verif_tetrad(tetrad,y1,Given_Metric)
        time_coord=Float64(y1[1])
        r=Float64(y1[3])
        dr=Float64(y1[4])
        phi=mod2pi(Float64(y1[5])) # longitude modulo 2*pi
        theta=Float64(y1[7])       # colatitude
        tau=Float64(y1[9])
        dtaudevpar_int=der1_tau_internal(y1,Given_Metric)
        # this is a quasi-Doppler for an ingoing spherical wave, with the Doppler counter on the S/C
        # can be outputted in a file, if needed
        Doppler=transmit_frequency_Doppler*(y1[2]-BF1) # in Hertz
        four_velocity_check=Float64(dtaudevpar_int-BF1)
        push!(four_velocity_check_plot,log10(abs(four_velocity_check)))
        push!(time_coord_plot,time_coord/86400.0)
        push!(r_plot,r)
        push!(dr_plot,dr)
        push!(Doppler_plot,Float64(Doppler)) 
        # X, Y, Z Cartesian coordinates from Boyer–Lindquist polar coordinates,
        # for function Kerr_Metric_Polar with as=Jstar/Mstar/c_light
        Z=r*cos(theta)
        X=sqrt(r^2+as^2)*sin(theta)*cos(phi) # this is not the conversion valid in Euclidean space
        Y=sqrt(r^2+as^2)*sin(theta)*sin(phi) # this is not the conversion valid in Euclidean space
        push!(Z_Plot,Z/1000.0)
        push!(X_Plot,X/1000.0)
        push!(Y_Plot,Y/1000.0)
        push!(phi_plot,phi/pi*180.0)
        push!(theta_plot,theta/pi*180.0)
        push!(tau_plot,tau/86400.0)
        push!(dif_tetrad_plot,Float64(log10(abs(dif_tetrad))))
        global evpar_count=evpar_count+evpar_step
    end
end
tc0=time_coord_plot[1]
taustart=tau_plot[1]
dif_tc_tau_plot=Float64[]
for i=1:length(time_coord_plot)
    time_coord_plot[i]=time_coord_plot[i]-tc0
    tau_plot[i]=tau_plot[i]-taustart
    push!(dif_tc_tau_plot,time_coord_plot[i]-tau_plot[i])
end
#----------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=plot(tau_plot,r_plot,"g",linewidth=1.0)
title("radius distance v. proper time _$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("radius (km)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_radius$(timetag).png"
savefig(filefig)
#----------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=scatter(tau_plot,phi_plot,marker=".",linewidths=0.2)
title("longitude v. proper time _$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("longitude (deg.)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_longitude$(timetag).png"
savefig(filefig)
#-----------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=scatter(tau_plot,theta_plot,marker=".",linewidths=0.2)
title("colatitude v. proper time _$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("colatitude (deg.)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_colatitude$(timetag).png"
savefig(filefig)
#-------------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=plot(tau_plot,dr_plot,"g",linewidth=1.0)
title("radial velocity v. proper time _$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("radial velocity (km/s)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_dr$(timetag).png"
savefig(filefig)
#------------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=plot(tau_plot,Doppler_plot,"g",linewidth=1.0)
title("Doppler_ingoing_wave_$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("Doppler (Hertz)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_Doppler$(timetag).png"
savefig(filefig)
#------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=plot(tau_plot,dif_tc_tau_plot,"g",linewidth=1.0)
title("time Coord - Proper time v. Proper time_$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("time Coord - Proper time (sec.) ",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_diftctau$(timetag).png"
savefig(filefig)
#------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=scatter(theta_plot,r_plot,marker=".",linewidths=0.2)
title("orbit(r,theta)_$(timetag)")
xlabel("theta",size=12)
ylabel("radius",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_orbit(r,theta)$(timetag).png"
savefig(filefig)
#------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=scatter(phi_plot,r_plot,marker=".",linewidths=0.2)
title("orbit(r,phi)_$(timetag)")
xlabel("phi",size=12)
ylabel("radius",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_orbit(r,phi)$(timetag).png"
savefig(filefig)
#-------------------------------------------------------------------------
figure(figsize=(10.5,9.5))
#plot_id=scatter([0.],[0.],[0.],marker=".",linewidth=12.0,"red")
plot_id=plot3D(X_Plot,Y_Plot,Z_Plot,linewidth=1.0)
title("orbit3D_$(timetag)")
xlabel("X (1000 km)",size=12)
ylabel("Y (1000 km)",size=12)
zlabel("Z (1000 km)",size=12)
display(plot_id)
filefig="display_3Dorbit$(timetag).png"
savefig(filefig)
#-------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=plot(tau_plot,dif_tetrad_plot,"g",linewidth=1.0)
title("FW tetrad check _$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("log10(error)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_FWtetrad_error$(timetag).png"
savefig(filefig)
#--------------------------------------------------------------------------
figure(figsize=(9.5,5.0))
plot_id=plot(tau_plot,four_velocity_check_plot,"g",linewidth=1.0)
title("four_velocity_check_$(timetag)")
xlabel("Proper time in days",size=12)
ylabel("log10(error)",size=12)
xticks(size=12)  
yticks(size=12)
display(plot_id)
filefig="display_four_velocity_check$(timetag).png"
savefig(filefig)
#--------------------------------------------------------------------------
prt3("timetag=",timetag)
msg=Call_System("cmd /c echo computation finished")
prt3(msg)
close(runlogfile)

end
