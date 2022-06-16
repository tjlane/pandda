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
            pos = (tuple(pos) if pos is not None else None),
            )

    def draw(self, obj):
        exec(self.as_cmd(obj))

    def as_cmd(self, obj, state=0):
        arg_string = ''
        for k, v in self.arg_dict.items():
            if v is not None: 
                if isinstance(v, str):
                    v = '"'+v.strip(' "')+'"'
                # else:
                #     v = str(v)
                arg_string += ', %s=%s' % (k, v)
        return 'cmd.pseudoatom("%s", state=%s%s)' % (obj, state, arg_string)


class Shape(object):

    def draw(self, obj):
        # would only ever be called inside a pymol session
        # so use print rather than logging
        print('Loading {shape} to {obj} '.format(shape=self.SHAPE, obj=obj))
        exec(self.as_cmd(obj))

    def as_cmd(self, obj, state=0):
        return 'cmd.load_cgo(%s, "%s", state=%s)' % (self.as_string(), obj, state)

    def as_string(self):
        raise Exception('Not Defined')


class Cylinder(Shape):

    SHAPE = 'Cylinder'

    def __init__(self, start, end, colors=None, radius=None):
        if (colors is None):
            colors = [[1,0,0],[1,0,0]]
        if (radius is None):
            radius = 0.1
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

    def __init__(self, centre, radius=None, color=None):
        if (color is None):
            color = [1,1,1]
        if (radius is None):
            radius = 1
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


class Cone(Shape):

    SHAPE = 'Cone'

    def __init__(self, start, end, colors=None, radius=None):
        if (colors is None):
            colors = [[1,0,0],[1,0,0]]
        if (radius is None):
            radius = 0.1
        self.start = start
        self.end = end
        self.colors = colors
        self.radius = radius

    def as_string(self):
        s = "["
        s += "CONE," + " %.6g, %.6g, %.6g," % tuple(self.start)
        s += " %.6g, %.6g, %.6g," % tuple(self.end)
        s += " %.6g, 0.0," % self.radius
        s += " %.6g, %.6g, %.6g," % tuple(self.colors[0])
        s += " %.6g, %.6g, %.6g," % tuple(self.colors[1])
        s += " 1.0, 1.0"
        s += "]"
        return s


class Arrow(object):

    SHAPE = 'Arrow'

    def __init__(self, start, end, colors=None, radius=None, head_length_fraction=0.4, head_radius_fraction=1.5):
        if (colors is None):
            colors = [[1,0,0],[1,0,0]]
        if (radius is None):
            radius = 0.1
        self.start = start
        self.end = end
        self.colors = colors
        self.radius = radius
        self.head_length_fraction = head_length_fraction
        self.head_radius_fraction = head_radius_fraction

        assert self.head_length_fraction >= 0.
        assert self.head_length_fraction <= 1.
        assert self.head_radius_fraction > 0.

    def make_components(self):

        body_length_fraction = (1.0 - self.head_length_fraction)

        start = tuple(self.start)
        end = tuple(self.end)
        mid = tuple(
            start[i] + body_length_fraction * (end[i] - start[i])
            for i in range(3)
            )

        c_start = list(self.colors[0])
        c_end = list(self.colors[1])
        c_mid = [
            c_start[i] + body_length_fraction * (c_end[i] - c_start[i])
            for i in range(3)
            ]

        return [
            Cylinder(
                start = start,
                end = mid,
                radius = self.radius,
                colors = [c_start, c_mid],
                ),
            Cone(
                start = mid,
                end = end,
                radius = (self.head_radius_fraction * self.radius),
                colors = [c_mid, c_end],
                ),
        ]

    def as_cmd(self, obj, state=0):
        return 'cmd.load_cgo(%s, "%s", state=%s)' % (self.as_string(), obj, state)

    def as_string(self):
        cpts = self.make_components()
        s = '['
        for c in cpts:
            s += c.as_string().strip('[] ')
            s += ','
        s = s.rstrip(', ') + ']'
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
