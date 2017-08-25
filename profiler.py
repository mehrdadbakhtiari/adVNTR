import time

from settings import *


def time_usage(func):
    def wrapper(*args, **kwargs):
        beg_ts = time.time()
        retval = func(*args, **kwargs)
        end_ts = time.time()
        function_name = func.__name__
        logging.info("%s executed in %fs" % (function_name, end_ts - beg_ts))
        return retval
    return wrapper
