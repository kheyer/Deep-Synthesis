# Streamlit confirm button with caching
# Source:
# https://gist.github.com/treuille/bc4eacbb00bfc846b73eec2984869645
# Written by Adrien Treuille
# https://gist.github.com/treuille

import streamlit as st
import collections
import functools
import inspect
import textwrap
from session_id import *

def cache_on_button_press(label, show_spinner=True):
    """Function decorator to memoize function executions.

    Parameters
    ----------
    label : str
        The label for the button to display prior to running the cached funnction.
    show_spinner : bool
        Whether to show the spinner when evaluating the cached function.

    Example
    -------
    This show how you could write a username/password tester:

    >>> @cache_on_button_press('Authenticate')
    ... def authenticate(username, password):
    ...     return username == "buddha" and password == "s4msara"
    ...
    ... username = st.text_input('username')
    ... password = st.text_input('password')
    ...
    ... if authenticate(username, password):
    ...     st.success('Logged in.')
    ... else:
    ...     st.error('Incorrect username or password')
    """
    def function_decorator(func):
        @functools.wraps(func)
        def wrapped_func(*args, **kwargs):
            @fancy_cache(allow_output_mutation=True, show_spinner=show_spinner,
                         unique_to_session=True, ttl=3600)
            def get_cache_entry(func, args, kwargs):
                class ButtonCacheEntry:
                    def __init__(self):
                        self.evaluated = False
                        self.return_value = None
                    def evaluate(self):
                        self.evaluated = True
                        self.return_value = func(*args, **kwargs)
                return ButtonCacheEntry()
            cache_entry = get_cache_entry(func, args, kwargs)
            if not cache_entry.evaluated:
                if st.button(label):
                    cache_entry.evaluate()
                else:
                    raise st.ScriptRunner.StopException
            return cache_entry.return_value
        return wrapped_func
    return function_decorator