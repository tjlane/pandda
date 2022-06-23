import giant.logs as lg
logger = lg.getLogger(__name__)

from libtbx import adopt_init_args


class MultiProcessWrapper(object):

    get_timeout = 0.05
    gc_interval = 60
    exit_flag = "STOP"

    def __init__(self,
            function,
            ):
        adopt_init_args(self, locals())

    def __call__(self,
            task_queue,
            done_queue,
            shutdown_event,
            ):

        import gc, time, traceback
        from multiprocessing.queues import Empty

        # Collect gc every gc_interval
        next_gc_time = time.time() + self.gc_interval

        while not shutdown_event.is_set():

            # Run GC if required (run at top of loop)
            if time.time() > next_gc_time:
                gc.collect()
                next_gc_time = time.time() + self.gc_interval

            # Get item or continue loop
            try:
                item = task_queue.get(block=True, timeout=self.get_timeout)
            except Empty:
                continue

            # Exit command has been issued
            if (item == self.exit_flag):
                break

            # Iterate through the arguments (tuples of (i,args))
            done_list = []
            for i, arg_dict in item:

                try:
                    output = (i, self.function(**arg_dict))
                except Exception:
                    output = (i, traceback.format_exc())

                done_list.append(output)

            # Append to output list
            done_queue.put(done_list)

            # Delete just in case helps gc
            del done_list


class RunParallelWithProgressBarUnordered(object):

    def __init__(self,
            function,
            n_cpus,
            max_chunksize = 10,
            keep_processes_open = False,
            ):

        processes = None
        task_queue = None
        done_queue = None
        shutdown_event = None

        assert max_chunksize >= 1

        adopt_init_args(self, locals())

    def __del__(self):
        self.close_processes()

    def __call__(self,
            arg_dicts,
            ):

        import tqdm
        import numpy

        # Initialise processes
        if (not self.processes):
            self.initialise_processes()

        try:

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

            # Get the results as they're ready
            results = []
            for i in range(n_chunks):
                complete = self.done_queue.get()
                results.extend(complete)
                pbar.update(len(complete))
            pbar.close()

            # Check that queues are empty
            self.check_queues_empty()

            # Sort output results
            results = [r[1] for r in sorted(results, key=lambda x: x[0])]

        except (Exception, KeyboardInterrupt) as e:
            self.close_processes()
            raise
        finally:
            pbar.close()

        # Close processes unless keep_open
        if (self.keep_processes_open is False):
            self.close_processes()

        return results

    def initialise_processes(self):

        from multiprocessing import Process, Queue, Event

        # Create input/output queues
        self.task_queue = Queue()
        self.done_queue = Queue()
        self.shutdown_event = Event()

        # Redirect termination signal while processes are spawned
        import signal
        sigint = signal.signal(signal.SIGINT, signal.SIG_IGN) # <- Set the signal handler to ignore
        self.processes = []
        for i in range(self.n_cpus):
            p = Process(
                target = MultiProcessWrapper(self.function),
                args = (self.task_queue, self.done_queue, self.shutdown_event),
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

        self.task_queue.close()
        self.task_queue.join_thread()

        self.done_queue.close()
        self.done_queue.join_thread()

        self.processes = None
        self.task_queue = None
        self.done_queue = None
        self.shutdown_event = None

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
