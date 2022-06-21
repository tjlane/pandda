import giant.logs as lg
logger = lg.getLogger(__name__)

from giant.utils import (
    make_sorted_dict,
    show_dict,
    merge_dicts,
    pretty_format_list,
    )


class DataCollator(object):

    def __init__(self):

        self.data = {}

    def __call__(self, update_dict):

        merge_dicts(
            master_dict = self.data,
            merge_dict = update_dict,
            )

        return self.data

    def __str__(self):

        s = json.dumps(
            self.data,
            indent = 2,
            )

        return s

    def get_sorted(self, recursive=True):

        return make_sorted_dict(
            self.data,
            recursive = recursive,
            )