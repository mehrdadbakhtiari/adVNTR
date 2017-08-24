import time

from settings import *


def time_usage(func):
    def wrapper(*args, **kwargs):
        beg_ts = time.time()
        retval = func(*args, **kwargs)
        end_ts = time.time()
        function_name = func.__name__
        with open(PROFILING_OUTPUT, 'a') as profiling_out:
            profiling_out.write("%s executed in %fs\n" % (function_name, end_ts - beg_ts))
        return retval
    return wrapper
