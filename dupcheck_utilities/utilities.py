import time
import sys


def log(message):
    """ Log messages to standard error. """
    sys.stderr.write(time.ctime() + ' --- ' + message + "\n")
    sys.stderr.flush()


def gc_content(s):
    """
    Calculate %GC in a sequence
    :param s: sequence for which to calculate %GC
    :return: float
    """
    total = 0
    for i in s:
        if any([i == 'G', i == 'g', i == 'C', i == 'c']):
            total += 1
    return total/len(s)