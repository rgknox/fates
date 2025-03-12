import numpy as np

def GetCosZ(ZONE,JULIAN,LAT,LON,LHOUR):

    #SUBROUTINE COSZENITH (LON,LATD,LHOUR,ZONE,JULIAN,CZENITH)
    #
    #     The purpose is to calculate the following:
    #        1)  Day angle (GAMMA)
    #        2)  Solar DEClination
    #        3)  Equation of time
    #        4)  Local apparent time
    #        5)  Hour angle
    #        6)  Cosine of zenith angle
    #
    #     All equations come from "An Introduction to
    #     Solar Radition" By Muhammad Iqbal, 1983.
    #
    # Integers:
    #     ZONE,                ! time zone (1-24) GMT=12
    #     JULIAN               ! julian day
    # Real
    #     CZENITH,             ! cosine of zenith angle (radians)
    #     DECd,                ! solar declination (degrees)
    #     DEC,                 ! solar declination (radians)
    #     et,                  ! equation of time (minutes)
    #     GAMMA,               ! day angle (radians)
    #     LATime,              ! local apparent time
    #     LCORR,               ! longitudical correction
    #     LHOUR,               ! local standard time
    #     LON,                 ! local longitude (deg)
    #     LLAT,                ! local latitude in radians
    #     LATD ,               ! local latitude in degrees
    #     LS,                  ! standard longitude (deg)
    #     OMEGAD,              ! omega in degrees
    #     OMEGA ,              ! omega in radians
    #     PI,                  ! universal constant PI [-]
    #     ZENITH,              ! zenith angle(radians)
    #     ZEND                 ! zenith angle(degress)
    #
    #     Neither ZENITH nor ZEND are necessary for this program.
    #     I originally used them as checks, and left them here in
    #     case anyone else had a use for them.
    #
    #     1)  Day angle GAMMA (radians) page 3

    PI= 3.141592              # universal constant PI
    GAMMA=2.*PI*(JULIAN-1)/365.

    #     2) Solar declination (assumed constant for a 24 hour period)  page 7
    #     in radians
    #
    DEC=(0.006918-0.399912*np.cos(GAMMA)+0.070257*np.sin(GAMMA) \
         -0.006758*np.cos(2.*GAMMA)+0.000907*np.sin(2.*GAMMA) \
         -0.002697*np.cos(3.*GAMMA)+0.00148*np.sin(3.*GAMMA))
    DECd=DEC*(180./PI)
    #
    #     maximum error 0.0006 rad (<3'), leads to error of less than 1/2 degree
    #     in ZENITH angle
    #     3)  Equation of time  page 11

    et=(0.000075+0.001868*np.cos(GAMMA)-0.032077*np.sin(GAMMA) \
        -0.014615*np.cos(2*GAMMA)-0.04089*np.sin(2*GAMMA))*229.18
    #
    #     4) Local apparent time  page 13
    #
    #     LS     standard longitude (nearest 15 degree meridian)
    #     LON     local longitude
    #     LHOUR  local standard time
    #     LATIME local apparent time
    #     LCORR  longitudunal correction (minutes)
    #
    
    LS=((ZONE-1)*15)-180.
    LCORR=4.*(LS-LON)*(-1.)
    LATIME=LHOUR+LCORR/60.+et/60.
    if (LATIME<0.):
        LATIME=LATIME+24.
    if (LATIME>24.):
        LATIME=LATIME-24.

    #     5) Hour angle OMEGA  page 15
    #
    #     hour angle is zero at noon, postive in the morning
    #     It ranges from 180 to -180
    #
    #     OMEGAD is OMEGA in degrees, OMEGA is in radians
    #
    OMEGAD=(LATIME-12.)*(-15.)
    OMEGA=OMEGAD*PI/180.
    #
    #     6)  Zenith angle page 15
    #
    #     CZENITH cosine of zenith angle (radians)
    #     LAT=local latitude in degrees
    #     LLAT=local latitude in radians
    #
    LLAT=LAT*PI/180.
    CZENITH=np.sin(DEC)*np.sin(LLAT)+np.cos(DEC)*np.cos(LLAT)*np.cos(OMEGA)
    CZENITH=np.max([0.,CZENITH])
    ZENITH=np.arcsin(CZENITH)
    ZEND=180.*ZENITH/PI
    
    return CZENITH
