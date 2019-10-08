# https://gist.github.com/tvst/6ef6287b2f3363265d51531c62a84f51
# Streamlit session specific caching
# Written by Thiago Teixeira
# https://gist.github.com/tvst

import streamlit.ReportThread as ReportThread
from streamlit.server.Server import Server


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