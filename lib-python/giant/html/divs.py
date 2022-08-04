class Div(object):

    __slots__ = (
        'id',
        'contents',
    )

    type = 'div'

    def __init__(self, **kw_args):

        # Ensure has contents attribute
        self.contents = []

        # Apply supplied args
        self.set(**kw_args)

    def append(self, other):
        self.contents.append(other)
        return other

    def extend(self, other):
        self.contents.extend(other)

    def set(self, **kw_args):
        for k, v in list(kw_args.items()):
            setattr(self, k, v)

    def get(self, key):
        assert key != 'contents'
        return getattr(self, key)

    def __getitem__(self, attr):
        return getattr(self, attr)

    def __setitem__(self, attr, value):
        setattr(self, attr, value)

    def update(self, dictionary):
        self.set(**dictionary)


class Block(Div):

    __slots__ = (
        'title',
        'text',
        'html',
        'image',
        'table',
        'footnote',
        'width',
        'colour',
        'classes',
        'styles',
        'max_width',
        'fancy_title',
        'title_size',
    )

    type = 'block'

    def __init__(
        self,
        width = 12,
        fancy_title = False,
        **kw_args
        ):

        # MUST BE FIRST
        super(Block, self).__init__(**kw_args)

        # Not captured by kw_args
        self.set(
            width = width,
            fancy_title = fancy_title,
        )
        

class Alert(Block):

    type = 'alert'


class ScrollX(Block):

    type = 'scroll'


class Panel(Block):

    __slots__ = Block.__slots__ + (
        'show',
    )

    type = 'panel'


class TabSet(Div):

    __slots__ = Div.__slots__ + (
        'title',
        'width',
        'title_size',
        'classes',
        'colour',
    )

    type = 'tabs'

    def set_active(
        self,
        tab = None,
        i_tab = 0,
        ):

        if not self.contents:
            return

        if tab is None:
            tab = self.contents[i_tab]

        assert tab in self.contents

        tab.active = True


class Tab(Block):

    __slots__ = Block.__slots__ + (
        'alt_title',
        'active',
    )

    type = 'tab'

    def __init__(
        self,
        alt_title = None,
        fancy_title = True,
        **kw_args
        ):

        super(Tab, self).__init__(**kw_args)

        if alt_title is None:
            alt_title = kw_args['title']

        self.set(
            alt_title = alt_title,
            fancy_title = fancy_title,
        )

