try:
    from pymol import cmd
    from pymol.cgo import *
except:
    pass

import collections


class Pseudoatom(object):


    def __init__(self, 
        selection = None,
        name = None,
        resn = None, 
        resi = None,
        chain = None,
        segi = None,
        elem = None,
        vdw = None,
        hetatm = None,
        b = None,
        q = None,
        color = None,
        label = None,
        pos = None,
        ):

        self.arg_dict = collections.OrderedDict(
            selection = selection,
            name = name,
            resn = resn, 
            resi = resi,
            chain = chain,
            segi = segi,
            elem = elem,
            vdw = vdw,
            hetatm = hetatm,
            b = b,
            q = q,
            color = color,
            label = label,
            pos = pos,
            )

    def draw(self, obj):
        exec(self.as_cmd(obj))

    def as_cmd(self, obj, state=0):
        arg_string = ''
        for k, v in self.arg_dict.iteritems():
            if v is not None: 
                if isinstance(v, str):
                    v = '"'+v.strip(' "')+'"'
                arg_string += ', %s=%s' % (k, v)
        return 'cmd.pseudoatom("%s", state=%s%s)' % (obj, state, arg_string)


class Shape(object):


    def draw(self, obj):
        print 'Loading %s to %s ' % (self.SHAPE, obj)
        exec(self.as_cmd(obj))

    def as_cmd(self, obj, state=0):
        return 'cmd.load_cgo(%s, "%s", state=%s)' % (self.as_string(), obj, state)

    def as_string(self):
        raise Exception('Not Defined')


class Cylinder(Shape):


    SHAPE = 'Cylinder'

    def __init__(self, start, end, colors=None, radius=None):
        if (colors is None):    colors = [[1,0,0]]*2
        if (radius is None):     radius  = 0.1
        self.start = start
        self.end = end
        self.colors = colors
        self.radius = radius

    def as_string(self):
        s = "["
        s += "CYLINDER," + " %.6g, %.6g, %.6g," % tuple(self.start)
        s += " %.6g, %.6g, %.6g," % tuple(self.end)
        s += " %.6g," % self.radius
        s += " %.6g, %.6g, %.6g," % tuple(self.colors[0])
        s += " %.6g, %.6g, %.6g" % tuple(self.colors[1])
        s += "]"
        return s


class Sphere(Shape):
    

    SHAPE = 'Sphere'

    def __init__(self, centre, radius=None, color=[1,1,1]):
        if (radius is None):    radius  = 1
        self.centre = centre
        self.radius = radius
        self.color = color

    def as_string(self):
        s = "["
        s += "COLOR,"
        s += " %.6g, %.6g, %.6g," % tuple(self.color)
        s += "SPHERE,"
        s += " %.6g, %.6g, %.6g," % tuple(self.centre)
        s += " %.6g" % self.radius
        s += "]"
        return s


class Alpha(Shape):


    SHAPE = 'None'

    def __init__(self, alpha_value=0.5):
        self.alpha = alpha_value

    def as_string(self):
        s = "["
        s += "ALPHA, %.6g" % self.alpha
        s += "]"
        return s
