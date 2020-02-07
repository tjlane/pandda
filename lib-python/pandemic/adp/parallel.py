import time, tqdm, signal
from libtbx import adopt_init_args

import traceback, multiprocessing, multiprocessing.pool


class MultiProcessWrapper:


    def __init__(self,
            function,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            task_queue,
            done_queue,
            ):

        for i, arg_dict in iter(task_queue.get, 'STOP'):

            try:
                output = (i, self.function(**arg_dict))
            except:
                tr = traceback.format_exc()
                if hasattr(self.function, 'log'):
                    self.function.log(tr)
                output = (i, tr)

            done_queue.put(output)


class RunParallelWithProgressBarUnordered:


    def __init__(self,
            function,
            n_cpus,
            chunksize = 10,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            arg_dicts,
            ):

        from multiprocessing import Process, Queue

        task_queue = Queue()
        done_queue = Queue()

        function = MultiProcessWrapper(self.function)

        # Redirect termination signal while processes are spawned
        sigint = signal.signal(signal.SIGINT, signal.SIG_IGN) # set the signal handler to ignore
        processes = []
        for i in range(self.n_cpus):
            p = Process(target=function, args=(task_queue, done_queue))
            p.start()
            processes.append(p)
        signal.signal(signal.SIGINT, sigint) # Replace the original signal handler

        for i, task in enumerate(arg_dicts):
            task_queue.put((i, task))

        pbar = tqdm.tqdm(total=len(arg_dicts), ncols=100)

        try:
            results = []
            for i in range(len(arg_dicts)):
                complete = done_queue.get()
                pbar.update(1)
                results.append(complete)
            pbar.close()
        except KeyboardInterrupt:
            for p in processes:
                p.terminate()
            raise
        finally:
            pbar.close()
            for p in processes:
                task_queue.put('STOP')

        # Sort output results
        results = [r[1] for r in sorted(results, key=lambda x: x[0])]

        return results
