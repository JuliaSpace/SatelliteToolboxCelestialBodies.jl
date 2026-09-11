## Description #############################################################################
#
# Constants used to compute the position and velocity of the Moon.
#
## References ##############################################################################
#
# [1] Meeus, J. (1998). Astronomical Algorithms. 2nd ed. Willmann-Bell, Inc, Richmond, VA.
#
############################################################################################

# The following tables were obtained from [1, pp. 339-341]. They are stored as tuples of
# static vectors so that the loops over the terms can be unrolled at compile time.

"""
    const _TAB_47A::NTuple{60, SVector{6, Int32}}

Periodic terms for the longitude (Σl) and distance (Σr) of the Moon, obtained from the
table 47.A in **[1, pp. 339-340]**.

Each element contains, in order, the multipliers of the fundamental arguments `D`, `M`,
`M´`, and `F` [-], the coefficient of the sine term for the longitude [10⁻⁶ °], and the
coefficient of the cosine term for the distance [m].

# References

- **[1]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.
"""
const _TAB_47A = (
#! format: off
    #                  D   M   M´  F        Σl         Σr
    #                 ─────────────── ────────  ─────────
    SVector{6, Int32}( 0,  0,  1,  0,  6288774, -20905355),
    SVector{6, Int32}( 2,  0, -1,  0,  1274027,  -3699111),
    SVector{6, Int32}( 2,  0,  0,  0,   658314,  -2955968),
    SVector{6, Int32}( 0,  0,  2,  0,   213618,   -569925),
    SVector{6, Int32}( 0,  1,  0,  0,  -185116,     48888),
    SVector{6, Int32}( 0,  0,  0,  2,  -114332,     -3149),
    SVector{6, Int32}( 2,  0, -2,  0,    58793,    246158),
    SVector{6, Int32}( 2, -1, -1,  0,    57066,   -152138),
    SVector{6, Int32}( 2,  0,  1,  0,    53322,   -170733),
    SVector{6, Int32}( 2, -1,  0,  0,    45758,   -204586),
    SVector{6, Int32}( 0,  1, -1,  0,   -40923,   -129620),
    SVector{6, Int32}( 1,  0,  0,  0,   -34720,    108743),
    SVector{6, Int32}( 0,  1,  1,  0,   -30383,    104755),
    SVector{6, Int32}( 2,  0,  0, -2,    15327,     10321),
    SVector{6, Int32}( 0,  0,  1,  2,   -12528,         0),
    SVector{6, Int32}( 0,  0,  1, -2,    10980,     79661),
    SVector{6, Int32}( 4,  0, -1,  0,    10675,    -34782),
    SVector{6, Int32}( 0,  0,  3,  0,    10034,    -23210),
    SVector{6, Int32}( 4,  0, -2,  0,     8548,    -21636),
    SVector{6, Int32}( 2,  1, -1,  0,    -7888,     24208),
    SVector{6, Int32}( 2,  1,  0,  0,    -6766,     30824),
    SVector{6, Int32}( 1,  0, -1,  0,    -5163,     -8379),
    SVector{6, Int32}( 1,  1,  0,  0,     4987,    -16675),
    SVector{6, Int32}( 2, -1,  1,  0,     4036,    -12831),
    SVector{6, Int32}( 2,  0,  2,  0,     3994,    -10445),
    SVector{6, Int32}( 4,  0,  0,  0,     3861,    -11650),
    SVector{6, Int32}( 2,  0, -3,  0,     3665,     14403),
    SVector{6, Int32}( 0,  1, -2,  0,    -2689,     -7003),
    SVector{6, Int32}( 2,  0, -1,  2,    -2602,         0),
    SVector{6, Int32}( 2, -1, -2,  0,     2390,     10056),
    SVector{6, Int32}( 1,  0,  1,  0,    -2348,      6322),
    SVector{6, Int32}( 2, -2,  0,  0,     2236,     -9884),
    SVector{6, Int32}( 0,  1,  2,  0,    -2120,      5751),
    SVector{6, Int32}( 0,  2,  0,  0,    -2069,         0),
    SVector{6, Int32}( 2, -2, -1,  0,     2048,     -4950),
    SVector{6, Int32}( 2,  0,  1, -2,    -1773,      4130),
    SVector{6, Int32}( 2,  0,  0,  2,    -1595,         0),
    SVector{6, Int32}( 4, -1, -1,  0,     1215,     -3958),
    SVector{6, Int32}( 0,  0,  2,  2,    -1110,         0),
    SVector{6, Int32}( 3,  0, -1,  0,     -892,      3258),
    SVector{6, Int32}( 2,  1,  1,  0,     -810,      2616),
    SVector{6, Int32}( 4, -1, -2,  0,      759,     -1897),
    SVector{6, Int32}( 0,  2, -1,  0,     -713,     -2117),
    SVector{6, Int32}( 2,  2, -1,  0,     -700,      2354),
    SVector{6, Int32}( 2,  1, -2,  0,      691,         0),
    SVector{6, Int32}( 2, -1,  0, -2,      596,         0),
    SVector{6, Int32}( 4,  0,  1,  0,      549,     -1423),
    SVector{6, Int32}( 0,  0,  4,  0,      537,     -1117),
    SVector{6, Int32}( 4, -1,  0,  0,      520,     -1571),
    SVector{6, Int32}( 1,  0, -2,  0,     -487,     -1739),
    SVector{6, Int32}( 2,  1,  0, -2,     -399,         0),
    SVector{6, Int32}( 0,  0,  2, -2,     -381,     -4421),
    SVector{6, Int32}( 1,  1,  1,  0,      351,         0),
    SVector{6, Int32}( 3,  0, -2,  0,     -340,         0),
    SVector{6, Int32}( 4,  0, -3,  0,      330,         0),
    SVector{6, Int32}( 2, -1,  2,  0,      327,         0),
    SVector{6, Int32}( 0,  2,  1,  0,     -323,      1165),
    SVector{6, Int32}( 1,  1, -1,  0,      299,         0),
    SVector{6, Int32}( 2,  0,  3,  0,      294,         0),
    SVector{6, Int32}( 2,  0, -1, -2,        0,      8752),
#! format: on
)

"""
    const _TAB_47B::NTuple{60, SVector{5, Int32}}

Periodic terms for the latitude (Σb) of the Moon, obtained from the table 47.B in
**[1, p. 341]**.

Each element contains, in order, the multipliers of the fundamental arguments `D`, `M`,
`M´`, and `F` [-], and the coefficient of the sine term for the latitude [10⁻⁶ °].

# References

- **[1]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.
"""
const _TAB_47B = (
#! format: off
    #                  D   M   M´  F        Σb
    #                 ─────────────── ────────
    SVector{5, Int32}( 0,  0,  0,  1,  5128122),
    SVector{5, Int32}( 0,  0,  1,  1,   280602),
    SVector{5, Int32}( 0,  0,  1, -1,   277693),
    SVector{5, Int32}( 2,  0,  0, -1,   173237),
    SVector{5, Int32}( 2,  0, -1,  1,    55413),
    SVector{5, Int32}( 2,  0, -1, -1,    46271),
    SVector{5, Int32}( 2,  0,  0,  1,    32573),
    SVector{5, Int32}( 0,  0,  2,  1,    17198),
    SVector{5, Int32}( 2,  0,  1, -1,     9266),
    SVector{5, Int32}( 0,  0,  2, -1,     8822),
    SVector{5, Int32}( 2, -1,  0, -1,     8216),
    SVector{5, Int32}( 2,  0, -2, -1,     4324),
    SVector{5, Int32}( 2,  0,  1,  1,     4200),
    SVector{5, Int32}( 2,  1,  0, -1,    -3359),
    SVector{5, Int32}( 2, -1, -1,  1,     2463),
    SVector{5, Int32}( 2, -1,  0,  1,     2211),
    SVector{5, Int32}( 2, -1, -1, -1,     2065),
    SVector{5, Int32}( 0,  1, -1, -1,    -1870),
    SVector{5, Int32}( 4,  0, -1, -1,     1828),
    SVector{5, Int32}( 0,  1,  0,  1,    -1794),
    SVector{5, Int32}( 0,  0,  0,  3,    -1749),
    SVector{5, Int32}( 0,  1, -1,  1,    -1565),
    SVector{5, Int32}( 1,  0,  0,  1,    -1491),
    SVector{5, Int32}( 0,  1,  1,  1,    -1475),
    SVector{5, Int32}( 0,  1,  1, -1,    -1410),
    SVector{5, Int32}( 0,  1,  0, -1,    -1344),
    SVector{5, Int32}( 1,  0,  0, -1,    -1335),
    SVector{5, Int32}( 0,  0,  3,  1,     1107),
    SVector{5, Int32}( 4,  0,  0, -1,     1021),
    SVector{5, Int32}( 4,  0, -1,  1,      833),
    SVector{5, Int32}( 0,  0,  1, -3,      777),
    SVector{5, Int32}( 4,  0, -2,  1,      671),
    SVector{5, Int32}( 2,  0,  0, -3,      607),
    SVector{5, Int32}( 2,  0,  2, -1,      596),
    SVector{5, Int32}( 2, -1,  1, -1,      491),
    SVector{5, Int32}( 2,  0, -2,  1,     -451),
    SVector{5, Int32}( 0,  0,  3, -1,      439),
    SVector{5, Int32}( 2,  0,  2,  1,      422),
    SVector{5, Int32}( 2,  0, -3, -1,      421),
    SVector{5, Int32}( 2,  1, -1,  1,     -366),
    SVector{5, Int32}( 2,  1,  0,  1,     -351),
    SVector{5, Int32}( 4,  0,  0,  1,      331),
    SVector{5, Int32}( 2, -1,  1,  1,      315),
    SVector{5, Int32}( 2, -2,  0, -1,      302),
    SVector{5, Int32}( 0,  0,  1,  3,     -283),
    SVector{5, Int32}( 2,  1,  1, -1,     -229),
    SVector{5, Int32}( 1,  1,  0, -1,      223),
    SVector{5, Int32}( 1,  1,  0,  1,      223),
    SVector{5, Int32}( 0,  1, -2, -1,     -220),
    SVector{5, Int32}( 2,  1, -1, -1,     -220),
    SVector{5, Int32}( 1,  0,  1,  1,     -185),
    SVector{5, Int32}( 2, -1, -2, -1,      181),
    SVector{5, Int32}( 0,  1,  2,  1,     -177),
    SVector{5, Int32}( 4,  0, -2, -1,      176),
    SVector{5, Int32}( 4, -1, -1, -1,      166),
    SVector{5, Int32}( 1,  0,  1, -1,     -164),
    SVector{5, Int32}( 4,  0,  1, -1,      132),
    SVector{5, Int32}( 1,  0, -1, -1,     -119),
    SVector{5, Int32}( 4, -1,  0, -1,      115),
    SVector{5, Int32}( 2, -2,  0,  1,      107),
#! format: on
)
