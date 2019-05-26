import time, tqdm, signal
from libtbx import adopt_init_args

import traceback, multiprocessing, multiprocessing.pool

class NonDaemonicProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class NonDaemonicPool(multiprocessing.pool.Pool):
    Process = NonDaemonicProcess


class MultiProcessWrapper:


    def __init__(self,
            function,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            arg_dict,
            ):
        i = arg_dict.pop('sort_value')
        try:
            return (i, self.function(**arg_dict))
        except:
            tr = traceback.format_exc()
            if hasattr(self.function, 'log'):
                self.function.log(tr)
                #self.function.log.log_file().write(clear_data=True)
            return (i, tr)


class RunParallelWithProgressBarUnordered:


    def __init__(self,
            n_cpus,
            chunksize = 10,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            function,
            arg_list,
            ):

        # Redirect termination signal
        sigint = signal.signal(signal.SIGINT, signal.SIG_IGN) # set the signal handler to ignore
        pool = NonDaemonicPool(self.n_cpus)
        signal.signal(signal.SIGINT, sigint) # Replace the original signal handler
        results = []
        cur_chunksize = min(self.chunksize, 1+int(float(len(arg_list)-1)/float(self.n_cpus)))
        try:
            pbar = tqdm.tqdm(total=len(arg_list), ncols=100)
            for complete in pool.imap_unordered(func=function, iterable=arg_list, chunksize=cur_chunksize):
                pbar.update(1)
                results.append(complete)
            pbar.close()
        except KeyboardInterrupt:
            pbar.close()
            pool.terminate()
            pool.join()
            raise
        # Close pool
        pool.close()
        pool.join()

        return results
