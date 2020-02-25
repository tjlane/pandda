import time, tqdm, signal
import numpy
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

        # Get `chunksize` objects from the queue
        for arg_dicts in iter(task_queue.get, 'STOP'):

            done_list = []

            # Iterate through the arguments (tuples of (i,args))
            for i, arg_dict in arg_dicts:

                try:
                    output = (i, self.function(**arg_dict))
                except:
                    output = (i, traceback.format_exc())
                done_list.append(output)

            done_queue.put(done_list)


class RunParallelWithProgressBarUnordered:


    def __init__(self,
            function,
            n_cpus,
            max_chunksize = 10,
            ):
        assert max_chunksize >= 1
        adopt_init_args(self, locals())

    def __call__(self,
            arg_dicts,
            ):

        n_tasks = len(arg_dicts)

        # Calculate actual chunksize (to efficiently use multi-CPUs)
        cpu_chunksize = numpy.ceil( n_tasks / self.n_cpus )
        chunksize = int(max(1, min(self.max_chunksize, cpu_chunksize)))

        from multiprocessing import Process, Queue

        task_queue = Queue()
        done_queue = Queue()

        function = MultiProcessWrapper(self.function)

        # Redirect termination signal while processes are spawned
        # v Set the signal handler to ignore v
        sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
        processes = []
        for i in range(self.n_cpus):
            p = Process(target=function, args=(task_queue, done_queue))
            p.start()
            processes.append(p)
        # v Replace the original signal handler v
        signal.signal(signal.SIGINT, sigint)

        # Number the tasks for sorting later and split into chunks
        task_list = list(enumerate(arg_dicts))
        n_chunks = 0
        for i in range(0, n_tasks, chunksize):
            chunk = task_list[i:i+chunksize]
            task_queue.put(chunk)
            n_chunks += 1

        pbar = tqdm.tqdm(total=n_tasks, ncols=100)

        try:
            results = []
            for i in range(n_chunks):
                complete = done_queue.get()
                pbar.update(len(complete))
                results.extend(complete)
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
