import giant.logs as lg
logger = lg.getLogger(__name__)

import os
from libtbx import adopt_init_args
from libtbx.utils import Sorry, Failure
from bamboo.html import png2base64src_maybe

from . import divs
import pandemic.resources

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


# Base class for html outputs
class HtmlSummary:

    # Permanent counter than links all subclasses!
    # Ensures all divs are uniquely numbered
    counter = Counter()

    debug = False
    embed_images = True

    @staticmethod
    def wrap_string(string, tag='p'):
        return '<'+tag+'>'+str(string)+'</'+tag+'>'

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
                    l = ':<strong>'.join(l.split(':')) + "</strong>"
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
    def image(cls, path):
        if (path is None) or not os.path.exists(path):
            path = pandemic.resources.NO_IMAGE_PATH
        if cls.embed_images is True:
            return png2base64src_maybe(path, print_on_missing=cls.debug)
        else:
            return path


class HtmlSummaryCollator(HtmlSummary):


    def __init__(self,
        title,
        alt_title = None,
        summaries = [],
        ):
        adopt_init_args(self, locals())

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


class HtmlSummaryConcatenator:


    def __init__(self,
            title,
            alt_title = None,
            summaries = [],
            ):
        adopt_init_args(self, locals())

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


class WriteHtmlSummaryTask:


    output_filename = 'results.html'

    def __init__(self,
        output_directory,
        verbose = False,
        ):
        output_filename_path = os.path.join(output_directory, self.output_filename)
        adopt_init_args(self, locals())

    def run(self,
        html_objects,
        ):

        # Things for the first page
        short_summaries = []
        main_summaries = []
        json_plots = []

        # Compile output lists
        for obj in html_objects:
            if hasattr(obj, 'short_summary'):
                short_summaries.extend(obj.short_summary())
            if hasattr(obj, 'main_summary'):
                main_summaries.extend(obj.main_summary())
            if hasattr(obj, 'json_plots'):
                json_plots.extend(obj.json_plots())

        self.make_output(
            overview_objects = short_summaries,
            tab_objects = main_summaries,
            json_plots = json_plots,
            )

        return self

    def make_output(self,
        overview_objects,
        tab_objects,
        json_plots = [],
        ):

        logger.subheading('Writing output HTML')

        # Get the html template
        from pandemic.html import PANDEMIC_HTML_ENV
        template = PANDEMIC_HTML_ENV.get_template('adp_summary.html')

        # Hard-coded for now...
        header_title = 'PanDEMIC ADP Summary'
        body_header = self.get_body_header()

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {
            'header_title' : header_title,
            'body_header' : body_header,
            'contents' : [],
            }
        # Jsons
        if json_plots:
            output_data['json_plots'] = json_plots
        # ===========================================================>

        # Create overview tab
        tab_set = divs.TabSet()
        tab_set.append(self.create_overview_tab(objects=overview_objects))
        tab_set.set_active()

        # Create other tabs
        tab_set.extend(tab_objects)

        output_data['contents'].append(tab_set)

        # Write out and format
        with open(self.output_filename_path, 'w') as out_html:
            out_html.write(template.render(output_data).encode( "utf-8" ))

        logger('Output HTML written to {}'.format(self.output_filename_path))

    def get_body_header(self):
        d = {
            'title' : 'Hierarchical Disorder Characterisation Summary',
            'sub_title' : 'Please cite: << Paper Title Here >>',
            'description' : 'N. Pearce & P. Gros, 2020.',
            'image' : HtmlSummary.image(pandemic.resources.FULL_LOGO_PATH),
        }

        return d

    def create_overview_tab(self, objects):

        tab = divs.Tab(
            id = 'overview',
            title = 'Hierarchical Disorder Parameterisation Overview',
            alt_title = 'Overview',
            contents = objects,
        )

        return tab

