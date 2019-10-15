# https://gist.github.com/tvst/6ef6287b2f3363265d51531c62a84f51
# Streamlit session specific caching
# Written by Thiago Teixeira
# https://gist.github.com/tvst

import streamlit as st
import streamlit.ReportThread as ReportThread
from streamlit.server.Server import Server
import time
import functools
import random
import string


def get_session_id():
    # Hack to get the session object from Streamlit.

    ctx = ReportThread.get_report_ctx()

    session = None
    session_infos = Server.get_current()._session_infos.values()

    for session_info in session_infos:
        if session_info.session._main_dg == ctx.main_dg:
            session = session_info.session

    if session is None:
        raise RuntimeError(
            "Oh noes. Couldn't get your Streamlit Session object"
            'Are you doing something fancy with threads?')

    return id(session)


# https://gist.github.com/treuille/f988f78c4610c78322d089eb77f74598
# Streamlit session specific cache decorator with expiration
# Written by Adrien Treuille
# https://gist.github.com/treuille

def fancy_cache(func=None, ttl=None, unique_to_session=False, **cache_kwargs):
    """A fancier cache decorator which allows items to expire after a certain time
    as well as promises the cache values are unique to each session.
    Parameters
    ----------
    func : Callable
        If not None, the function to be cached.
    ttl : Optional[int]
        If not None, specifies the maximum number of seconds that this item will
        remain in the cache.
    unique_to_session : boolean
        If so, then hash values are unique to that session. Otherwise, use the default
        behavior which is to make the cache global across sessions.
    **cache_kwargs
        You can pass any other arguments which you might to @st.cache
    """
    # Support passing the params via function decorator, e.g.
    # @fancy_cache(ttl=10)
    if func is None:
        return lambda f: fancy_cache(
            func=f,
            ttl=ttl,
            unique_to_session=unique_to_session,
            **cache_kwargs
        )
    
    # This will behave like func by adds two dummy variables.
    dummy_func = st.cache(
        func = lambda ttl_token, session_token, *func_args, **func_kwargs: \
            func(*func_args, **func_kwargs),
        **cache_kwargs)
    
    # This will behave like func but with fancy caching.
    @functools.wraps(func)
    def fancy_cached_func(*func_args, **func_kwargs):
        # Create a token which changes every ttl seconds.
        ttl_token = None
        if ttl is not None:
            ttl_token = int(time.time() / ttl)
        
        # Create a token which is unique to each session.
        session_token = None
        if unique_to_session:
            session_token = get_session_id()
        
        # Call the dummy func
        return dummy_func(ttl_token, session_token, *func_args, **func_kwargs)
    return fancy_cached_func