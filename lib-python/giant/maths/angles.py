import numpy

def rad2deg(a):
    return a*180.0/numpy.pi
def deg2rad(a):
    return a*numpy.pi/180.0

def angle_difference(a1, a2, deg=True, abs_val=False):
    """Calculate the difference between two angles"""

    if deg is False:
        a1 = rad2deg(a1)
        a2 = rad2deg(a2)

    d = (a2-a1+180.0)%360.0-180.0

    if abs_val:
        d = numpy.abs(d)

    if deg is False:
        return deg2rad(d)
    else:
        return d
