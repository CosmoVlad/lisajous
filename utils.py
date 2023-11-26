import numpy as np

import sympy as sp
sp.init_printing()


## conversion between two pairs of angles characterizing the orbital momentum of a source: 
## w.r.t. the Sun (thetaL, phiL) and w.r.t. sky plane (iota -- inclination, Phi -- polarization angle)

def Ln(thetaS,phiS,thetaL,phiL):
    
    return sp.cos(thetaL)*sp.cos(thetaS) + sp.sin(thetaL)*sp.sin(thetaS)*sp.cos(phiL-phiS)
    
def convert_angle(angle):
    
    # convert the angle range [-pi,pi] (the output of atan2) to [0,2*pi]
    
    if not isinstance(angle, np.ndarray):
        if angle < 0:
            return angle + 2*np.pi
        return angle
    
    angle_copy = angle.copy()
    angle_copy[angle_copy < 0] += 2*np.pi
    return angle_copy

def source_polarization(thetaS, phiS, thetaL, phiL):
    
    return sp.atan2(
        sp.sin(thetaL)*sp.sin(phiL-phiS), 
        sp.sin(thetaL)*sp.cos(phiL-phiS)*sp.cos(thetaS) - sp.cos(thetaL)*sp.sin(thetaS)
    )

def source_inclination(thetaS, phiS, thetaL, phiL):
    
    return sp.acos(Ln(thetaS,phiS,thetaL,phiL))

def PhiL(thetaS, phiS, iota, Phi):
    
    return sp.atan2(
        sp.cos(iota)*sp.sin(thetaS)*sp.sin(phiS) + sp.sin(iota)*(
            sp.sin(Phi)*sp.cos(phiS) + sp.cos(Phi)*sp.sin(phiS)*sp.cos(thetaS)
        ),
        sp.cos(iota)*sp.sin(thetaS)*sp.cos(phiS) - sp.sin(iota)*(
            sp.sin(Phi)*sp.sin(phiS) - sp.cos(Phi)*sp.cos(phiS)*sp.cos(thetaS)
        )
    )

def ThetaL(thetaS, phiS, iota, Phi):
    
    return sp.acos(
        sp.cos(iota)*sp.cos(thetaS) - sp.sin(iota)*sp.sin(thetaS)*sp.cos(Phi)
    )
    
## LISA pattern functions in terms of (iota, Phi)
## time is measured in yr

def phi_t(t):
    
    return 2*sp.pi*t

def expr_cos(t, theta, phi):
    
    return sp.cos(theta)/2 - np.sqrt(3)/2 * sp.sin(theta)*sp.cos(phi_t(t) - phi)

def Lz(t, thetaL, phiL):
    
    return expr_cos(t, thetaL, phiL)

def expr_cos_thetaS(t, thetaS, phiS):
    
    return expr_cos(t, thetaS, phiS)


def expr_phiS(t, thetaS, phiS):
    
    return phi_t(t) + sp.atan((np.sqrt(3)*sp.cos(thetaS) + sp.sin(thetaS)*sp.cos(phi_t(t)-phiS))\
                             /(2*sp.sin(thetaS)*sp.sin(phi_t(t)-phiS)))

def nLz(t, thetaS, phiS, thetaL, phiL):
    
    A = sp.cos(thetaL)*sp.sin(thetaS)*sp.sin(phiS) - sp.cos(thetaS)*sp.sin(thetaL)*sp.sin(phiL)
    B = sp.cos(thetaS)*sp.sin(thetaL)*sp.cos(phiL) - sp.cos(thetaL)*sp.sin(thetaS)*sp.cos(phiS)
    
    return sp.sin(thetaL)*sp.sin(thetaS)*sp.sin(phiL-phiS)/2\
                - np.sqrt(3)/2 * sp.cos(phi_t(t)) * A\
                - np.sqrt(3)/2 * sp.sin(phi_t(t)) * B


def polarization_angle(t, thetaS, phiS, Phi):
    
    z1 = sp.sqrt(3)/2 * sp.cos(thetaS) * sp.cos(phi_t(t)-phiS) + sp.Rational(1,2) * sp.sin(thetaS)
    z2 = sp.sqrt(3)/2 * sp.sin(phi_t(t)-phiS)
    
    return sp.atan2((-z1*sp.cos(Phi) - z2*sp.sin(Phi)),(z1*sp.sin(Phi) - z2*sp.cos(Phi)))

def Fplus(t, thetaS, phiS, iota, Phi, phase_shift=0):
    
    cosThetaS = expr_cos_thetaS(t, thetaS, phiS)
    pol_angle = polarization_angle(t, thetaS, phiS, Phi)
    PhiS = expr_phiS(t, thetaS, phiS)
    PhiS -= phase_shift
    
    return (1+cosThetaS**2)/2 * sp.cos(2*PhiS)*sp.cos(2*pol_angle)\
                        - cosThetaS * sp.sin(2*PhiS)*sp.sin(2*pol_angle)

def Fcross(t, thetaS, phiS, iota, Phi, phase_shift=0):
    
    cosThetaS = expr_cos_thetaS(t, thetaS, phiS)
    pol_angle = polarization_angle(t, thetaS, phiS, Phi)
    PhiS = expr_phiS(t, thetaS, phiS)
    PhiS -= phase_shift
    
    return (1+cosThetaS**2)/2 * sp.cos(2*PhiS)*sp.sin(2*pol_angle)\
                        + cosThetaS * sp.sin(2*PhiS)*sp.cos(2*pol_angle)


