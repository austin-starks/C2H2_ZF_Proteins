"""
A collection of helper functions for project.py.

This file is used to store helper functions for project.py. These functions are
used to edit filenames, display the time elapsed, as well as other functions.

INSTANCE ATTRIBUTES:
    last_second:        the last second to be displayed by time_elapsed [int]

Author: Austin Starks
Date: July 4, 2019
"""

import os.path
from time import time
last_second = -1


def time_elapsed(start_time):
    """
    Prints the amount of time that has elapsed from a given start time. If the
    message is not the empty string, it will display a message as well as the
    elapsed time.

    Parameter start_time: the time the timer begins at
    Preconditon: start_time is an int

    Parameter msg: the message to display
    Precondtion: msg is a string.
    """
    time_ = round(time() - start_time)
#    print(str(time_))
    seconds = time_ % 60
    minutes = time_ // 60
    if minutes < 10:
        str_minutes = "0" + str(minutes)
    else:
        str_minutes = str(minutes)
    if seconds < 10:
        str_seconds = "0" + str(seconds)
    else:
        str_seconds = str(seconds)
    global last_second
    if seconds != last_second:
        easy_print("Elapsed time", str_minutes +
                   ":" + str_seconds, same_line=False)
        # same_line is False to make printing on Atom console easier
        # same_line should be True if using console
        last_second = seconds


def edit_filename(filename, location, foldername=''):
    """
    Edits a filename to be correct.

    This function takes in a filename and gets the correct directory of the file.

    Parameter filename: the name of the file
    Preconditon: filename is a string

    Parameter location: the locaxtion of the file.
    Preconditon: location is a "child", "parent", "sibling", "same"

    Parameter foldername: the name of the sibling/child folder
    Preconditon: foldername is a string
    """
    assert type(filename) == str, "Filename must be a string"
    assert type(location) == str, "Location must be a string"
    assert location is "child" or location is "parent" or location is "sibling" or \
        location is "same"
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    if location == "child":
        filename = os.path.join(fileDir, foldername + '/' + filename)
    elif location == "parent":
        filename = os.path.join(fileDir, "../" + filename)
    elif location == "sibling":
        filename = os.path.join(fileDir, '../' + foldername + '/' + filename)
    elif location == "same":
        pass
    return filename


def easy_print(message, value, same_line=False):
    """
    Returns: An error message about the given value.

    Parameter message: The error message to display
    Precondition: message is a string

    Parameter value: The value that caused the error
    Precondition: NONE (value can be anything)

    Parameter same_line: says whether to print the message on the same
    line. Used for loops
    Preconditon: same_line is a bool
    """
    assert type(message) == str, "Message must be a string"
    assert type(same_line) == bool, "Same_line must be a bool"
    if not same_line:
        print(message + ': ' + str(value))
    else:
        print("\r" + message + ': ' + str(value), end="", flush=True)


class DoesNotExistError(Exception):
    """
    This Error class is thrown when a file does not exist.
    """
    pass


class IncorrectError(Exception):
    """
    This Error class is thrown when a nucleotide does not match the expected
    nucleotide.
    """
    pass
