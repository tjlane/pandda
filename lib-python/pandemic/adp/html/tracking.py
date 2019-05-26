import os
from libtbx import adopt_init_args
from pandemic.adp.html import HtmlSummary


class TrackingHtmlSummary(HtmlSummary):


    def __init__(self,
        tracking_object,
        ):
        adopt_init_args(self, locals())

    def short_summary(self):
        output = [
            {
                'type'       : 'panel',
                'title'      : 'Optimisation Summary',
                'contents'   : [
                    {
                        'width' : 6,
                        'image' : self.image(
                            self.tracking_object.tracking_png1,
                            ),
                        },
                    {
                        'width' : 6,
                        'image' : self.image(
                            self.tracking_object.tracking_png2,
                            ),
                        },
                    ],
                },
            ]
        return output


