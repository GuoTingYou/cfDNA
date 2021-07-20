import os
import sys
import time
from functools import wraps, partial

def timemain(func):
    @wraps(func)
    def decorator(*args, **kwargs):
        print(os.path.abspath(os.path.realpath(sys.modules[func.__module__].__file__)))
        print("\n{:-^81s}".format(
              " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
        ))

        func(*args, **kwargs) # main function here

        print("{:-^81s}\n".format(
              " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
        ))
        print("Practice_makes_perfect")
    return decorator

def timetask(func):
    @wraps(func)
    def decorator(*args, **kwargs):
        t0 = time.time()
        print("{0} [{1}]".format(time.strftime("%X", time.localtime()), func.__name__))
        result = func(*args, **kwargs)
        print("{0} [{1}]".format(time.strftime("%X", time.localtime()), func.__name__))
        t1 = time.time() - t0
        print("Finished [{1}] consumed: {0} Secs.".format(round(t1, 2), func.__name__))
        return result
    return decorator

def bounded_classmethod(meth=None, classes:list=None):
    if meth is None:
        return partial(bounded_classmethod, classes=classes)
    for cls in classes:
        setattr(cls, meth.__name__, meth)
    return meth

def argParse_template(func):
    @wraps(func)
    def decorator(N):
        parser = func(N)
        if len(sys.argv) <= N:
            parser.print_help()
            raise IndexError(f"Not enough arguments were given. Required {N}")
        return parser.parse_args()
    return decorator
