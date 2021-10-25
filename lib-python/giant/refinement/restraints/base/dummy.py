from .collection import (
    RestraintsCollection,
    )


class DummyRestraintMaker(object):

    name = "DummyRestraintMaker"

    def __init__(self, name=None, *args, **kwargs):

        if name is not None:
            self.name = str(name)

    def __str__(self):

        s_ = (
            'Task: {name}\n'
            '| This task does nothing.\n'
            '`---->'
            ).format(
                name = self.name,
                )

        return s_.strip()

    def __call__(self, hierarchy):
        return RestraintsCollection()
