import streamlit as st

from src.common import *


params = page_setup(page="main")


st.markdown(
    """
# ICAnnoLncRNA
"""
)

save_params(params)
