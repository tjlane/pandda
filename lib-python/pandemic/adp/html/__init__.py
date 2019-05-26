import os
from libtbx import adopt_init_args
from bamboo.html import png2base64src_maybe


def as_html_summary_maybe(obj):
    if hasattr(obj, 'as_html_summary'):
        return obj.as_html_summary()
    else:
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
    def format_summary(text, width=12, type='alert', colour='info'):
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
            p_dict = {'type'   : type,
                      'width'  : width,
                      'text'   : '<br>'.join(fmt_lines),
                      'colour' : colour}
            if title is not None:
                p_dict['title'] = title
            output.append(p_dict)
        return output

    def short_summary(self): 
        return []
    def main_summary(self): 
        return []
    def json_plots(self): 
        return []

    def image(self, path):
        if self.embed_images is True:
            return png2base64src_maybe(path, print_on_missing=self.debug)
        else:
            return path


class HtmlSummaryCollator(HtmlSummary):


    def __init__(self,
        title,
        alt_title = None,
        summaries = [],
        ):
        if alt_title is None: alt_title = title
        adopt_init_args(self, locals())

    def main_summary(self):

        content_list = []
        for s in self.summaries:
            content_list += s.main_summary()
        for t in content_list:
            if not t.get('alt_title'):
                t['alt_title'] = t['title']

        if True or (len(content_list) > 1):
            content_list[0]['active'] = True
            contents = [{'type' : 'tabs', 'contents' : content_list}]
        else:
            contents = content_list

        output = {
            'title' : self.title,
            'alt_title' : self.alt_title,
            'contents' : contents,
        }

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


class WriteHtmlSummaryTask:


    output_filename = 'results.html'

    def __init__(self,
        output_directory,
        verbose = False,
        log = None
        ):
        if log is None: log = Log()
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
            header = 'PanDEMIC ADP Output Summary',
            title = 'PanDEMIC ADP Parameterisation',
            introduction = 'Hierarchical ADP parameterisation summary from pandemic.adp',
            overview_objects = short_summaries,
            tab_objects = main_summaries,
            json_plots = json_plots,
            )

        return self

    def make_output(self,
        header,
        title,
        introduction,
        overview_objects,
        tab_objects,
        json_plots = [],
        ):

        log = self.log
        log.subheading('Writing output HTML')

        # Get the html template
        from pandemic.html import PANDEMIC_HTML_ENV
        template = PANDEMIC_HTML_ENV.get_template('adp_summary.html')

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {}
        output_data['header'] = header
        output_data['title'] = title
        output_data['introduction'] = introduction
        # ===========================================================>
        # Construct the tabs
        output_data['contents'] = []

        # Create overview tab
        tab_set = {'type':'tabs', 'contents':[]}
        tab_set['contents'].append(self.create_overview_tab(objects=overview_objects))
        tab_set['contents'][-1]['active'] = True

        # Create other tabs
        tab_set['contents'].extend(tab_objects)

        output_data['contents'].append(tab_set)

        # Jsons
        if json_plots:
            output_data['json_plots'] = json_plots

        # Write out and format
        with open(self.output_filename_path, 'w') as out_html:
            out_html.write(template.render(output_data).encode( "utf-8" ))

        log('Output HTML written to {}'.format(self.output_filename_path))


    def create_overview_tab(self, objects):

        tab = {
            'id'        : 'overview',
            'alt_title' : 'Overview',
            'title'     : 'Hierarchical Disorder Parameterisation Overview',
            'contents'  : objects,
            }

        return tab

