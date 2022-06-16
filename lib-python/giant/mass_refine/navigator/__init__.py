
from giant.mass_refine import (
    MassRefineFileSystem,
    )

def coot_customisation():
    set_nomenclature_errors_on_read("ignore")
    set_recentre_on_read_pdb(0)
    set_show_symmetry_master(1)
    set_symmetry_shift_search_size(2)
    set_colour_map_rotation_for_map(0.0)
    set_colour_map_rotation_on_read_pdb(0.0)
    set_colour_map_rotation_on_read_pdb_flag(1)
    set_colour_map_rotation_on_read_pdb_c_only_flag(1)
    #set-stop-scroll-iso-map 0

class GuiPart(object):

    def __init__(self):

        self.labels = {}
        self.buttons = {}
        self.objects = {}

    def __call__(self):

        bdc_adjuster = gtk.Adjustment(
            lower = 0.0,
            upper = 1.0,
            value = float(starting_value),
            page_size = 0.01,
            step_incr = 0.01,
            page_incr = 0.1,
            )




class DatasetTracker(object):

    def __init__(self, file_system):

        self.file_system = file_system

        self._build()

    def _build(self):

        self.dataset_meta_table = self.get_dataset_table()

    def get_dataset_table(self):

        dataset_dicts = {}

        for d_lab, d_fs in self.file_system.datasets.get_all_as_dict():

            d_info = d_fs.get_meta()

            dataset_dicts[d_lab] = d_info

        return pd.DataFrame.from_dict(
            dataset_dicts,
            orient = 'index',
            )


class MassRefineNavigator(object):

    def __init__(self, dir_path):

        self.mass_refine_fs = MassRefineFileSystem(
            work_folder = dir_path,
            )

        self.dataset_table = self.build_dataset_table(
            file_system = self.mass_refine_fs,
            )

        self.update()

    def __call__():

        self.gui = NavigatorGui()


if __name__ == '__main__':

    #############################################################################################
    #
    # CHANGE COOT SETTINGS
    #
    #############################################################################################

    coot_setup()

    try:
        coot_customisation()
    except:
        pass

    #############################################################################################
    #
    # RUN
    #
    #############################################################################################

    # args = parser.parse_args(sys.argv)
    args, unknown = parser.parse_known_args(sys.argv)

    working_directory = pl.Path(
        args.pandda_directory
        # os.getcwd()
        )

    try:

        splash = SplashScreen()

        from pandda.inspect.io import (
            GetPanddaInspectInputOutputFiles,
            )

        get_io_files = GetPanddaInspectInputOutputFiles(
            input_directory = (
                working_directory
                ),
            output_directory = (
                working_directory / 'inspect'
                ),
            mode = args.mode,
            )

        inspector = PanddaInspect(
            working_directory = working_directory,
            files_dict = get_io_files(),
            mode = args.mode,
            )

        gui = inspector()

        splash.show_menu()

        splash.window.connect("destroy", gui.launch)
        # splash.window.connect("destroy_event", gui.launch)

        # gui.launch()

    except MissingFile as e:

        modal_msg(
            'Missing files! Are you in the right directory?\n\nNo {}\n\npandda.inspect will now close'.format(
                str(e)
                )
            )

        sys.exit()

