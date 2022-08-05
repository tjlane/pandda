import giant.logs as lg
logger = lg.getLogger(__name__)

import os

from . import divs

div_class_hash = dict(
    none = divs.Block,
    block = divs.Block,
    alert = divs.Alert,
    panel = divs.Panel,
)

def as_html_summary_maybe(obj):
    if obj is None:
        return None
    if hasattr(obj, 'as_html_summary'):
        return obj.as_html_summary()
    return None

def as_html_summaries_maybe(tasks):
    out = []
    for t in tasks:
        s = as_html_summary_maybe(t)
        if s is not None:
            out.append(s)
    return out


class Counter(object):
    def __init__(self, start=0):
        self.i = start
    def next(self):
        self.i += 1
        return self.i


class ImageEmbedder(object):

    def __init__(self, embed=False, relative_to=None):

        self.embed = embed
        self.relative_to = relative_to

    def __call__(self, image_path):

        # Replace with NO_IMAGE if necessary
        if (image_path is None):
            import pandemic.resources
            image_path = pandemic.resources.NO_IMAGE_PATH_ADP
        elif not os.path.exists(image_path) and (self.embed is True):
            import pandemic.resources
            image_path = pandemic.resources.NO_IMAGE_PATH_ADP

        # Read image and return as string
        if (self.embed is True):
            from giant.html import png2base64src_maybe
            return png2base64src_maybe(
                image_path, 
                print_on_missing = False,
                )

        # Create relative path
        if self.relative_to is not None:
            image_path = os.path.relpath(
                path = image_path, 
                start = self.relative_to,
                )

        return image_path


# Base class for html outputs
class HtmlSummary(object):

    # Permanent counter than links all subclasses!
    # Ensures all divs are uniquely numbered
    counter = Counter()

    debug = False
    embed_images = True

    # Link files relative to this path
    output_dir = None

    @staticmethod
    def wrap_string(string, tag='p'):
        return '<{t}>{s}</{t}>'.format(
            t = tag,
            s = str(string),
            )

    @staticmethod
    def format_summary(text, width=12, type='alert', colour=None, classes=None):
        div_obj_class = div_class_hash[type]
        paragraphs = text.split('\n\n')
        output = []
        for p in paragraphs:
            lines = p.strip('\n ').split('\n')
            if lines[0].startswith('>'):
                title = lines.pop(0).strip('>\n ')
            else:
                title = None
            fmt_lines = []
            for l in lines:
                if l.count(':')==1:
                    l = ':<strong>'.join(l.split(':')) + '</strong>'
                l = l.replace('\t','&emsp;')
                fmt_lines.append(l)
            obj = div_obj_class(
                width = width,
                contents = [divs.Block(text=l) for l in fmt_lines],
            )
            if title is not None:
                obj.title = title
            if colour is not None:
                obj.colour = colour
            if classes is not None:
                obj.classes = classes
            output.append(obj)
        return output

    def short_summary(self):
        return []
    def main_summary(self):
        return []
    def json_plots(self):
        return []

    @classmethod
    def image(cls, path, embed=None):

        if embed is None: 
            embed = cls.embed_images

        # Replace with NO_IMAGE if necessary
        if (path is None) or not os.path.exists(path):
            import pandemic.resources
            path = pandemic.resources.NO_IMAGE_PATH_ADP

        # Read image and return as string
        if (embed is not False):
            from giant.html import png2base64src_maybe
            return png2base64src_maybe(path, print_on_missing=cls.debug)

        # Create relative path if contained in output folder
        if (cls.output_dir is not None) and path.startswith(cls.output_dir):
            path = os.path.relpath(path=path, start=cls.output_dir)

        return path


class HtmlSummaryCollator(HtmlSummary):


    def __init__(self,
        title,
        alt_title = None,
        summaries = [],
        ):

        self.title = title
        self.alt_title = alt_title
        self.summaries = summaries

    def append(self, summary):
        self.summaries.append(summary)

    def main_summary(self):

        content_list = []
        for s in self.summaries:
            content_list.extend(s.main_summary())
        for t in content_list:
            assert isinstance(t, divs.Tab)
            assert (t.alt_title is not None)

        contents = divs.TabSet(
            contents = content_list,
        )
        contents.set_active()

        output = divs.Tab(
            title = self.title,
            alt_title = self.alt_title,
            contents = [contents],
        )

        return [output]

    def short_summary(self):
        l = []
        for s in self.summaries:
            l += s.short_summary()
        return l

    def json_plots(self):
        l = []
        for s in self.summaries:
            l += s.json_plots()
        return l


class HtmlSummaryConcatenator(object):


    def __init__(self,
            title,
            alt_title = None,
            summaries = [],
            ):
        
        self.title = title
        self.alt_title = alt_title
        self.summaries = summaries

    def append(self, summary):
        self.summaries.append(summary)

    def main_summary(self):
        collated = divs.Tab(
            title = self.title,
            alt_title = self.alt_title,
        )
        for s in self.summaries:
            s_summary = s.main_summary()
            if len(s_summary) == 0:
                continue
            for other_summary in s_summary:
                collated.extend(other_summary.contents)
        if collated is None:
            return []
        return [collated]

    def short_summary(self):
        l = []
        for s in self.summaries:
            l += s.short_summary()
        return l

    def json_plots(self):
        l = []
        for s in self.summaries:
            l += s.json_plots()
        return l
