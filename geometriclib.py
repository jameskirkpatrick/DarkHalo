"""geometriclib.py 
 functions used to work out geometric stuff
"""
from math import sqrt, cos, sin,pi
from numpy import arctan2

""" modulus
returns the modulus of a 2D vector"""
def modulus(coord):
    return sqrt(coord[0]**2 + coord[1]**2)

"""angle 
gets the angle of a 2D vector wrt to the origin
"""
def angle(coord):
    return arctan2(coord[1],coord[0])    

"""diff
take the difference of two 2D tuples
"""
def diff(something, minussomething):
    return (something[0]-minussomething[0], something[1]-minussomething[1])
 
"""rotateEllipticity
determines the tangential component of the ellipticity to a certain angle
"""
def rotateEllipticity(e, angle):
    eTangential = -e[0] *cos(2.0*angle) - e[1] * sin(2.0*angle)
    eCross      = e[0] *sin(2.0*angle) + e[1] * cos(2.0*angle)
    return (eTangential,eCross)

    
"""putAngleInNormalRange
returns an angle in the range -pi to pi
"""
def putAngleInNormalRange(angle):
    if angle > pi :
        return angle - 2.*pi
    elif angle < -pi :
        return angle + 2.*pi
    else:
        return angle   