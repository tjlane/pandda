import tqdm
import numpy
from libtbx import adopt_init_args

import traceback

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
            keep_processes_open = False,
            ):

        processes = None
        task_queue = None
        done_queue = None

        assert max_chunksize >= 1

        adopt_init_args(self, locals())

    def __del__(self):
        self.close_processes()

    def __call__(self,
            arg_dicts,
            ):

        # Initialise processes
        if (not self.processes):
            self.initialise_processes()

        # Number of actual jobs
        n_tasks = len(arg_dicts)

        # Calculate used chunksize (to efficiently use multi-CPUs)
        cpu_chunksize = numpy.ceil( n_tasks / self.n_cpus )
        chunksize = int(max(1, min(self.max_chunksize, cpu_chunksize)))

        # Check that task queue and done queues are empty
        self.check_queues_empty()
        self.check_processes_alive()

        # Number the tasks for sorting later and split into chunks
        task_list = list(enumerate(arg_dicts)) # <- TASKS NUMBERED HERE
        n_chunks = 0
        for i in range(0, n_tasks, chunksize):
            chunk = task_list[i:i+chunksize]
            self.task_queue.put(chunk)
            n_chunks += 1

        # Progress bar
        pbar = tqdm.tqdm(total=n_tasks, ncols=100)

        try:
            # Get the results as they're ready
            results = []
            for i in range(n_chunks):
                complete = self.done_queue.get()
                results.extend(complete)
                pbar.update(len(complete))
            pbar.close()

            # Check that queues are empty
            self.check_queues_empty()

        except KeyboardInterrupt:
            self.close_processes()
            raise
        except Exception as e:
            self.close_processes()
            raise e
        finally:
            pbar.close()

        # Sort output results
        results = [r[1] for r in sorted(results, key=lambda x: x[0])]

        # Close processes unless keep_open
        if (self.keep_processes_open is False):
            self.close_processes()

        return results

    def initialise_processes(self):

        from multiprocessing import Process, Queue

        # Create input/output queues
        self.task_queue = Queue()
        self.done_queue = Queue()

        # Redirect termination signal while processes are spawned
        import signal
        sigint = signal.signal(signal.SIGINT, signal.SIG_IGN) # <- Set the signal handler to ignore
        self.processes = []
        for i in range(self.n_cpus):
            p = Process(
                target = MultiProcessWrapper(self.function),
                args = (self.task_queue, self.done_queue),
            )
            p.start()
            self.processes.append(p)
        signal.signal(signal.SIGINT, sigint) # <- Replace the original signal handler

    def close_processes(self):

        if (self.processes is None):
            return

        # Pass flag to stop the processes
        for p in self.processes:
            if p.is_alive():
                self.task_queue.put('STOP')

        # Join all processes
        for p in self.processes:
            p.join(timeout=60)

        # Terminate the processes if still alive
        for p in self.processes:
            if p.is_alive():
                p.terminate()

        self.processes = None
        self.task_queue = None
        self.done_queue = None

    def check_queues_empty(self):
        try:
            if self.task_queue.empty() is False:
                raise Failure('Task queue is not empty')
            if self.done_queue.empty() is False:
                raise Failure('Done queue is not empty')
        except Exception as e:
            self.close_processes()
            raise e

    def check_processes_alive(self):
        try:
            for p in self.processes:
                if not p.is_alive():
                    raise Failure('One or more processes has been terminated prematurely.')
        except Exception as e:
            self.close_processes()
            raise e
