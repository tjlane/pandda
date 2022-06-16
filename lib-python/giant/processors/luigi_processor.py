import json
import luigi

from .luigi_sge import (
    SGEJobTask
    )

from .base_processors import (
    BaseProcessor,
    )


# TODO make LuigiCondor ? 

def random_string(length=30):
    import random, string

    return ''.join(
        random.choice(string.ascii_letters) for _ in range(length)
        )


class LuigiSGETask(SGEJobTask):

    func = luigi.Parameter()
    output_path = luigi.Parameter()

    def work(self):
        self.func()

    def output(self):
        return luigi.LocalTarget(
            str(self.output_path)
            )


class LuigiFunctionWrapper(object):

    def __init__(self, 
        function, 
        shared_tmp_dir, # random hex path
        output_file_suffix = None,
        ):

        self.function = (
            function
            )

        self.output_filename = (
            str(
                random_string(length=30)
                ) + str(
                '' if output_file_suffix is None else output_file_suffix
                )
            )

        self.output_path = (
            pl.Path(shared_tmp_dir) / output_filename 
            )

    def __call__(self):

        output = self.function()

        if (output is None) or (self.output_file_suffix is None):
            output_str = 'done'
        elif self.output_path.suffix == '.json':
            output_str = json.dumps(
                output
                )
        # pickle?
        else:
            output_str = 'done'

        with open(output_path) as fh:

            fh.write(output_str)

        return None


class ProcessorLuigiSGE(BaseProcessor):

    def __init__(self,
        n_nodes, 
        n_cpus_per_node,
        h_vmem, # default NONE do not rename v_mem_per_node
        m_mem_free, # default NONE do not rename free_mem_per_cpu?
        shared_tmp_dir,
        parallel_env="smp", # 
        output_file_suffix = None,
        ):

        self.n_nodes = n_nodes
        self.n_cpus_per_node = n_cpus_per_node

        self.h_vmem = h_vmem
        self.m_mem_free = m_mem_free

        self.shared_tmp_dir = shared_tmp_dir
        self.parallel_env = parallel_env

    def __call__(self,
        funcs,
        result_loaders = None,
        ):

        wrappers = [
            LuigiWrapper(f)
            for f in funcs
            ]

        tasks = [
            Task(
                func = w,
                output_path = w.output_path, # well-defined as using wrappers
                shared_tmp_dir = self.shared_tmp_dir, # "/dls/science/groups/i04-1/conor_dev/pandda/lib-python/pandda/pandda_analyse_dask/luigi_test",
                parallel_env = self.parallel_env,
                n_cpu = self.n_cpus_per_node,
                run_locally = False, # debug!
                h_vmem = self.h_vmem,
                m_mem_free = self.m_mem_free,
                )
            for w in wrappers 
            ]

        luigi.build(
            tasks,
            local_scheduler = True,
            workers = self.n_nodes,
            detailed_summary = False,
            )

        if result_loaders:
            results = [
                r()
                for r in result_loaders
                ]
        else:
            results = [
                None 
                for w in wrappers
                ]

        return results


