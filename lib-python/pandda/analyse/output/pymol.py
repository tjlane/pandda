import giant.logs as lg
logger = lg.getLogger(__name__)

import os
from giant.dispatcher import Dispatcher


class MakePymolOutputImages:

    output_front = "all_events_front.png"
    output_back = "all_events_back.png"
    output_script = "all_events_pymol.py"

    def __init__(self, 
        output_directory,
        run_script=True,
        ):

        self.output_directory = output_directory

        self.output_files = {
            'events_front' : str(
                self.output_directory / self.output_front
                ),
            'events_back' : str(
                self.output_directory / self.output_back
                ),
            }
            
        self.run_script = run_script

    def __call__(self, 
        reference_structure,
        event_dicts,
        site_dicts,
        ):

        if (reference_structure is None):
            raise ValueError('must supply reference_structure')

        output_script = self.make_pymol_script(
            reference_structure = reference_structure,
            event_dicts = event_dicts,
            site_dicts = site_dicts,
            )

        if self.run_script is True:
            try: 
                self.run_pymol_script(
                    script = str(output_script),
                    )
            except Exception as e: 
                logger(str(e))
                raise

        return self.output_files

    def make_pymol_script(self, 
        reference_structure, 
        event_dicts,
        site_dicts,
        ):

        site_counts = {
            s['site_num'] : 0 
            for s in site_dicts
        }
        for e in event_dicts:
            site_counts[e['site_num']] += 1

        ###

        pymol_str =  '# Mark the identified sites on the protein\n'
        pymol_str += 'from pymol import cmd\n'
        pymol_str += 'from pymol.cgo import *\n'
        pymol_str += 'cmd.load("{}", "reference")\n'.format(
            os.path.relpath(
                reference_structure, 
                start = str(self.output_directory),
                )
            )
        pymol_str += 'cmd.show_as("cartoon", "reference")\n'
        pymol_str += 'cmd.color("cyan", "reference")\n'

        for event in event_dicts:
            lab = 'event'
            com = tuple(event['xyz_centroid_ref'])
            pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=0.5)\n'.format(lab, com)
            pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
            pymol_str += 'cmd.color("blue", "{}")\n'.format(lab)

        for site in site_dicts:
            site_num = site['site_num']
            # Only print the site if it has more than one event
            if site_counts[site_num] == 1:
                continue
            lab = 'site_{}'.format(site_num)
            com = tuple(site['xyz_centroid'])
            pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=2.5)\n'.format(lab, com)
            pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
            pymol_str += 'cmd.label("{}", "{}")\n'.format(lab, site_num)
            pymol_str += 'cmd.color("deepteal", "{}")\n'.format(lab)
            pymol_str += 'cmd.set("label_color", "white", "{}")\n'.format(lab)

        # Set label things...
        pymol_str += 'cmd.set("label_size", 25)\n'
        pymol_str += 'cmd.set("label_position", (0,0,4))\n'
        pymol_str += 'cmd.bg_color(color="white")\n'

        # Write as python script
        pymol_script = str(
            self.output_directory / self.output_script
            )

        with open(pymol_script, 'w') as fh:
            fh.write(pymol_str)

        return pymol_script

    def run_pymol_script(self, script):

        pymol_str =  '# Load the protein representation and output images of sites\n'
        pymol_str += 'run {}\n'.format(
            os.path.relpath(
                script,
                start = str(self.output_directory),
                )
            )
        pymol_str += 'set ray_opaque_background, off\n'
        pymol_str += 'set specular, off\n'
        pymol_str += 'orient\n'
        pymol_str += 'png {}, width=1200, dpi=300, ray=0\n'.format(
            os.path.relpath(
                self.output_files['events_front'], 
                start = str(self.output_directory)
                )
            )
        pymol_str += 'turn y, 180\n'
        pymol_str += 'png {}, width=1200, dpi=300, ray=0\n'.format(
            os.path.relpath(
                self.output_files['events_back'], 
                start = str(self.output_directory)
                )
            )
        pymol_str += 'quit'

        script_pml = (
            str(script) + '.pml'
            )

        with open(script_pml, 'w') as fh:
            fh.write(pymol_str)

        # Change into directory as script runs off of relative paths
        cur_dir = os.getcwd()
        os.chdir(str(self.output_directory))

        try:
            c = Dispatcher('pymol')
        except ValueError as e:
            logger('Pymol is not available: {}'.format(str(e)))
            return

        c.extend_args([
            '-k', '-q', '-c', os.path.relpath(script_pml, start=str(self.output_directory))
            ])

        logger.subheading('Generating pymol images')
        logger(str(c))

        c.run()

        os.chdir(cur_dir)

        return