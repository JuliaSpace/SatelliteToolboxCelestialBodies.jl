# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Compute the Sun position and velocity.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] The Astronomical Almanac for the year 2000 (p. C24).
#
#   [2] http://aa.usno.navy.mil/faq/docs/SunApprox.php
#
#   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#       Microcosm Press, Hawthorne, CA.
#
#   [4] The Astronomical Almanac for the year 2006.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sun_position_mod, sun_velocity_mod, sun_position_el

"""
    sun_position_mod(jd_tdb::Number) -> SVector{3, Float64}
    sun_position_mod(date_tdb::DateTime) -> SVector{3, Float64}

Compute the Sun position represented in the IAU-76/FK5 MOD (mean-equator, mean-equinox of
date) at the Julian Day `jd_tdb` or `date_tdb`. The input time must be represented in the
Barycentric Dynamical Time (TDB). The algorithm was adapted from [1, pp. 277-279].

!!! Note
    This function performs all the computations using `Float64` due to the necessary
    precision.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
sun_position_mod(date_tdb::DateTime) = sun_position_mod(datetime2julian(date_tdb))

function sun_position_mod(jd_tdb::Number)
    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34t_tdb

    # Convert Ms to [rad] and limit to the interval [0,2π].
    Ms = mod2pi(deg2rad(Ms))

    # Compute auxiliary variables.
    sin_Ms,  cos_Ms  = sincos(Ms)

    # Let's compute `sin(2Ms)` and `cos(2Ms)` more efficiently.
    sin_2Ms = 2 * sin_Ms * cos_Ms
    cos_2Ms = cos_Ms * cos_Ms - sin_Ms * sin_Ms

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic latitude of the Sun [deg].
    λ_e = λ_m + 1.914_666_471sin_Ms + 0.019_994_643sin_2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2t_tdb

    # Convert λ_e and ϵ to [rad] and limit to the interval [0,2π].
    λ_e = mod2pi(deg2rad(λ_e))
    ϵ   = mod2pi(deg2rad(ϵ))

    # Auxiliary variables.
    sin_ϵ   , cos_ϵ   = sincos(ϵ)
    sin_λ_e , cos_λ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = (1.000_140_612 - 0.016_708_617cos_Ms - 0.000_139_589cos_2Ms) * ASTRONOMICAL_UNIT

    # Compute the Sun vector represented in the Mean Equinox of Date (MOD).
    s_mod = @SVector [r * cos_λ_e, r * cos_ϵ * sin_λ_e, r * sin_ϵ * sin_λ_e]

    return s_mod
end

"""
    sun_velocity_mod(jd_tdb::Number) -> SVector{3, Float64}
    sun_velocity_mod(date_tdb::DateTime) -> SVector{3, Float64}

Compute the Sun velocity measured and represented in the IAU-76/FK5 MOD (mean-equator,
mean-equinox of date) at the Julian Day `jd_tdb` or `date_tdb`. The input time must be
represented in the Barycentric Dynamical Time (TDB). The algorithm was obtained by computing
the time derivative of the Sun position in [1, p. 277-279].

!!! note
    This function performs all the computations using `Float64` due to the necessary
    precision.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
sun_velocity_mod(date_tdb::DateTime) = sun_velocity_mod(datetime2julian(date_tdb))

function sun_velocity_mod(jd_tdb::Number)
    # Convert centuries to seconds.
    cen2s = 1 / (36525.0 * 86400.0)

    # Sun position
    # ======================================================================================

    # TODO: Check if we can split this algorithm in a new function.
    #
    # This is the same algorithm as in `sun_position_mod`. However, we must have access to
    # some internal variables. Hence, we copied the entire function for now.

    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34t_tdb

    # Convert Ms to [rad] and limit to the interval [0,2π].
    Ms = mod2pi(deg2rad(Ms))

    # Compute auxiliary variables.
    sin_Ms,  cos_Ms  = sincos(Ms)

    # Let's compute `sin(2Ms)` and `cos(2Ms)` more efficiently.
    sin_2Ms = 2 * sin_Ms * cos_Ms
    cos_2Ms = cos_Ms * cos_Ms - sin_Ms * sin_Ms

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic latitude of the Sun [deg].
    λ_e = λ_m + 1.914_666_471sin_Ms + 0.019_994_643sin_2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2t_tdb

    # Convert λ_e and ϵ to [rad] and limit to the interval [0,2π].
    λ_e = mod2pi(deg2rad(λ_e))
    ϵ   = mod2pi(deg2rad(ϵ))

    # Auxiliary variables.
    sin_ϵ   , cos_ϵ   = sincos(ϵ)
    sin_λ_e , cos_λ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = (1.000_140_612 - 0.016_708_617cos_Ms - 0.000_139_589cos_2Ms) * ASTRONOMICAL_UNIT

    # Sun velocity
    # ==========================================================================

    # Compute the required time derivatives [rad/s].
    ∂λ_m = deg2rad(36_000.771)   * cen2s
    ∂Ms  = deg2rad(35_999.05034) * cen2s
    ∂ϵ   = deg2rad(-0.0130042)  * cen2s
    ∂r   = (+0.016708617sin_Ms * ∂Ms + 2 * 0.000139589cos_2Ms * ∂Ms) * ASTRONOMICAL_UNIT
    ∂λ_e = ∂λ_m + deg2rad(1.914666471cos_Ms * ∂Ms + 2 * 0.019994643cos_2Ms * ∂Ms)

    # Compute the Sun velocity vector represented in the Mean Equinox of Date
    # (MOD).
    vsun_mod = @SVector [
        ∂r * cos_λ_e -                                    r * sin_λ_e * ∂λ_e,
        ∂r * cos_ϵ * sin_λ_e - r * sin_ϵ * sin_λ_e * ∂ϵ + r * cos_ϵ * cos_λ_e * ∂λ_e,
        ∂r * sin_ϵ * sin_λ_e + r * cos_ϵ * sin_λ_e * ∂ϵ + r * sin_ϵ * cos_λ_e * ∂λ_e
    ]

    return vsun_mod
end

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Functions to Compute Sun Position In Local Frame

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                                  Aux. Functions
################################################################################

# Experimenting with accuracy of results by varying the method of computation of `t`
function give_t(yt, mt, Day, UT)
    opt = 1
    if opt == 1
        t = ((365.25*(yt-2000)) + (30.6001*(mt+1)) - (0.01*yt) + Day) + 0.0416667*UT - 21958.0;
    elseif opt == 2
        t = convert(Float64, floor(365.25*(yt-2000)) + floor(30.6001*(mt+1)) - floor(0.01*(yt)) + Day)
            + 0.0416667*UT - 21958.0;
    end
    return t
end

# Implementation of Algorithm1
function sunpos_alg1(
    UT::AbstractFloat,
    Day::Integer,
    Month::Integer,
    Year::Integer,
    Dt::AbstractFloat,
    Longitude::AbstractFloat,
    Latitude::AbstractFloat,
    Pressure::AbstractFloat,
    Temperature::AbstractFloat,
)

    (Month <= 2) ? (mt = Month + 12; yt = Year - 1) : (mt = Month; yt = Year)

    # Algorithm1
    t = give_t(yt, mt, Day, UT)
    te = t + 1.1574e-5*Dt;
    wte = 0.017202786*te;

    s1 = sin(wte);
    c1 = cos(wte);
    s2 = 2.0*s1*c1;
    c2 = (c1+s1)*(c1-s1);

    RightAscension = -1.38880 + 1.72027920e-2*te + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2
        + 1.525e-2*c2;
    RightAscension = rem(RightAscension, 2*π);
    Declination = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2;
    HourAngle = 1.75283 + 6.3003881*t + Longitude - RightAscension;
    HourAngle = rem(HourAngle + π, 2*π) - π;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    (ep > 0.0) ? (De = 0.08422*Pressure/((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)))) : (De = 0.0)

    Zenith = π/2 - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    return (RightAscension, Declination, HourAngle, Zenith, Azimuth)
end

# Implementation of Algorithm2
function sunpos_alg2(
    UT::AbstractFloat,
    Day::Integer,
    Month::Integer,
    Year::Integer,
    Dt::AbstractFloat,
    Longitude::AbstractFloat,
    Latitude::AbstractFloat,
    Pressure::AbstractFloat,
    Temperature::AbstractFloat,
)

    (Month <= 2) ? (mt = Month + 12; yt = Year - 1) : (mt = Month; yt = Year)

    # Algorithm2
    t = give_t(yt, mt, Day, UT)
    te = t + 1.1574e-5*Dt;

    wte = 0.017202786*te;

    s1 = sin(wte);
    c1 = cos(wte);
    s2 = 2.0*s1*c1;
    c2 = (c1+s1)*(c1-s1);
    s3 = s2*c1 + c2*s1;
    c3 = c2*c1 - s2*s1;
    s4 = 2.0*s2*c2;
    c4 = (c2+s2)*(c2-s2);

    RightAscension = -1.38880 + 1.72027920e-2*te + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2
        + 1.525e-2*c2 + 1.33e-3*s3 + 3.8e-4*c3 + 7.3e-4*s4 + 6.2e-4*c4;
    RightAscension = rem(RightAscension, 2*π);
    Declination = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2 + 1.50e-3*s3
        - 2.58e-3*c3 + 6e-5*s4 - 1.3e-4*c4;
    HourAngle = 1.75283 + 6.3003881*t + Longitude - RightAscension;
    HourAngle = rem(HourAngle + π, 2*π) - π;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    (ep > 0.0) ? (De = 0.08422*Pressure/((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)))) : (De = 0.0)

    Zenith = π/2 - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    return (RightAscension, Declination, HourAngle, Zenith, Azimuth)
end

# Implementation of Algorithm3
function sunpos_alg3(
    UT::AbstractFloat,
    Day::Integer,
    Month::Integer,
    Year::Integer,
    Dt::AbstractFloat,
    Longitude::AbstractFloat,
    Latitude::AbstractFloat,
    Pressure::AbstractFloat,
    Temperature::AbstractFloat,
)

    (Month <= 2) ? (mt = Month + 12; yt = Year - 1) : (mt = Month; yt = Year)

    # Algorithm3
    t = give_t(yt, mt, Day, UT)
    te = t + 1.1574e-5*Dt;

    wte = 0.0172019715*te;

    λ = -1.388803 + 1.720279216e-2*te + 3.3366e-2*sin(wte - 0.06172)
        + 3.53e-4*sin(2.0*wte - 0.1163);
    ϵ = 4.089567e-1 - 6.19e-9*te;

    sl = sin(λ);
    cl = cos(λ);
    se = sin(ϵ);
    ce = sqrt(1-se*se);

    RightAscension = atan(sl*ce, cl);
    (RightAscension < 0.0) ? (RightAscension += 2*π) : nothing
    Declination = asin(sl*se);
    HourAngle = 1.7528311 + 6.300388099*t + Longitude - RightAscension;
    HourAngle = rem(HourAngle + π, 2*π) - π;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    (ep > 0.0) ? (De = 0.08422*Pressure/((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)))) : (De = 0.0)

    Zenith = π/2 - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    return (RightAscension, Declination, HourAngle, Zenith, Azimuth)
end

# Implementation of Algorithm4
function sunpos_alg4(
    UT::AbstractFloat,
    Day::Integer,
    Month::Integer,
    Year::Integer,
    Dt::AbstractFloat,
    Longitude::AbstractFloat,
    Latitude::AbstractFloat,
    Pressure::AbstractFloat,
    Temperature::AbstractFloat,
)

    (Month <= 2) ? (mt = Month + 12; yt = Year - 1) : (mt = Month; yt = Year)

    # Algorithm4
    t = give_t(yt, mt, Day, UT)
    te = t + 1.1574e-5*Dt;

    wte = 0.0172019715*te;

    L = 1.752790 + 1.720279216e-2*te + 3.3366e-2*sin(wte - 0.06172)
        + 3.53e-4*sin(2.0*wte - 0.1163);

    ν = 9.282e-4*te - 0.8;
    Dlam = 8.34e-5*sin(ν);
    λ = L + π + Dlam;

    ϵ = 4.089567e-1 - 6.19e-9*te + 4.46e-5*cos(ν);

    sl = sin(λ);
    cl = cos(λ);
    se = sin(ϵ);
    ce = sqrt(1-se*se);

    RightAscension = atan(sl*ce, cl);
    (RightAscension < 0.0) ? (RightAscension += 2*π) : nothing
    Declination = asin(sl*se);
    HourAngle = 1.7528311 + 6.300388099*t + Longitude - RightAscension + 0.92*Dlam;
    HourAngle = rem(HourAngle + π, 2*π) - π;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    (ep > 0.0) ? (De = 0.08422*Pressure/((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)))) : (De = 0.0)

    Zenith = π/2 - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    return (RightAscension, Declination, HourAngle, Zenith, Azimuth)
end

# Implementation of Algorithm5
function sunpos_alg5(
    UT::AbstractFloat,
    Day::Integer,
    Month::Integer,
    Year::Integer,
    Dt::AbstractFloat,
    Longitude::AbstractFloat,
    Latitude::AbstractFloat,
    Pressure::AbstractFloat,
    Temperature::AbstractFloat,
)

    (Month <= 2) ? (mt = Month + 12; yt = Year - 1) : (mt = Month; yt = Year)

    # Algorithm5
    t = give_t(yt, mt, Day, UT)
    te = t + 1.1574e-5*Dt;

    wte = 0.0172019715*te;

    s1 = sin(wte);
    c1 = cos(wte);
    s2 = 2.0*s1*c1;
    c2 = (c1+s1)*(c1-s1);
    s3 = s2*c1 + c2*s1;
    c3 = c2*c1 - s2*s1;

    L = 1.7527901 + 1.7202792159e-2*te + 3.33024e-2*s1 - 2.0582e-3*c1 + 3.512e-4*s2
        - 4.07e-5*c2 + 5.2e-6*s3 - 9e-7*c3 - 8.23e-5*s1*sin(2.92e-5*te)
        + 1.27e-5*sin(1.49e-3*te - 2.337) + 1.21e-5*sin(4.31e-3*te + 3.065)
        + 2.33e-5*sin(1.076e-2*te - 1.533) + 3.49e-5*sin(1.575e-2*te - 2.358)
        + 2.67e-5*sin(2.152e-2*te + 0.074) + 1.28e-5*sin(3.152e-2*te + 1.547)
        + 3.14e-5*sin(2.1277e-1*te - 0.488);

    ν = 9.282e-4*te - 0.8;
    Dlam = 8.34e-5*sin(ν);
    λ = L + π + Dlam;

    ϵ = 4.089567e-1 - 6.19e-9*te + 4.46e-5*cos(ν);

    sl = sin(λ);
    cl = cos(λ);
    se = sin(ϵ);
    ce = sqrt(1-se*se);

    RightAscension = atan(sl*ce, cl);
    (RightAscension < 0.0) ? (RightAscension += 2*π) : nothing
    Declination = asin(sl*se);
    HourAngle = 1.7528311 + 6.300388099*t + Longitude - RightAscension + 0.92*Dlam;
    HourAngle = rem(HourAngle + π, 2*π) - π;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    (ep > 0.0) ? (De = 0.08422*Pressure/((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)))) : (De = 0.0)

    Zenith = π/2 - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    return (RightAscension, Declination, HourAngle, Zenith, Azimuth)
end

# Main Function
"""
    sun_position_el(
        JD::Number,
        Longitude::Real=0.0,
        Latitude::Real=0.0,
        Pressure::Real=1.0,
        Temperature::Real=20.0,
        flag::Char='l',
        algorithm::Char='5'
    )

Compute the Sun position represented in the Local Horizon Reference System at the
Julian Day `JD`, Longitude `Longitude`, Latitude `Latitude`, Atmospheric pressure `Pressure`,
Ambient Temperature `Temperature`, Output Flag `flag` and Algorithm `algorithm`.

Inputs:

JD: Julian Day;
Longitude and Latitude of the observer, in degrees, WSG84;
Pressure in atm;
Temperature in Celsius degrees.


Outputs:

Equatorial system: flag -e
    RightAscension, in radians;
    Declination, in radians;

Local Coordinates: flag -l
    HourAngle, in radians;
    Zenith, in radians;
    Azimuth, in radians;

flag -a: all outputs

"""
function sun_position_el(
    JD::Real,
    Longitude::Real=0.0,
    Latitude::Real=0.0,
    Pressure::Real=1.0,
    Temperature::Real=20.0,
    flag::Char='l',
    algorithm::Char='5',
)
    # Get time data from Julian Date `JD`
    Year, Month, Day, h, m, s = jd_to_date(JD)
    UT = h + m/60.0 + s/3600.0
    Dt = get_Δat(JD) + 32.184

    # Convert from degrees to radians
    Longitude = Longitude*π/180
    Latitude = Latitude*π/180

    # Switch algorithm based on the algorithm flag `algorithm`
    if algorithm == '1'
        Posn = sunpos_alg1(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature)
    elseif algorithm == '2'
        Posn = sunpos_alg2(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature)
    elseif algorithm == '3'
        Posn = sunpos_alg3(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature)
    elseif algorithm == '4'
        Posn = sunpos_alg4(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature)
    elseif algorithm == '5'
        Posn = sunpos_alg5(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature)
    end

    # Switch output based on the output flag `flag`
    if flag == 'e'
        # Equtorial System
        return Posn[1:2]
    elseif flag == 'l'
        # Local Coordinates
        return Posn[3:5]
    elseif flag == 'a'
        # All
        return Posn
    end
end
