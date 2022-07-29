
import copy

from . import divs


class DivFormatter(object):

    def __init__(self,
        div_class = divs.Block,
        *args, **kwargs
        ):

        self.div_class = div_class
        self.args = args
        self.kwargs = kwargs

    def __call__(self):

        return self.div_class()


class CitationFormatter(DivFormatter):


    def __call__(self,
        title, 
        journal,
        year,
        authors,
        link = None,
        div_class = divs.Block,
        ):

        text_lines = [
            self.format_title(title),
            self.format_authors(authors),
            self.format_journal_year(journal, year),
            ]
            
        if link is not None:
            text_lines.append(
                self.format_link(
                    link = link,
                    link_text = link,
                    )
                )

        args, kwargs = self.get_args()

        main_div = div_class(
            text = '<br>'.join(text_lines),
            *args, **kwargs
            )

        return main_div

    def get_args(self):

        args = copy.deepcopy(self.args)
        kwargs = copy.deepcopy(self.kwargs)

        return args, kwargs

    def format_title(self, title):

        return '<strong>{}</strong>'.format(title)

    def format_journal_year(self, journal, year):

        journal_year = '{journal} ({year})'.format(
            journal = journal,
            year = year,
            )

        return '<small><strong>{}</strong></small>'.format(journal_year)

    def format_authors(self, authors):

        if isinstance(authors, list):
            if len(authors) == 1:
                authors = authors[0]
            else: 
                authors = (
                    ', '.join(authors[:-1]) + 
                    ' and ' + 
                    authors[-1]
                    )
        else: 
            authors = str(authors)

        return '<small>{}</small>'.format(authors)

    def format_link(self, link, link_text):

        link_text = (
            '<a href="{link}"><small>{link_text}</small></a>'.format(
                link = link,
                link_text = link_text,
                )
            )

        return link_text

