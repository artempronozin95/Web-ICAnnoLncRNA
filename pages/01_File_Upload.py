from pathlib import Path

import streamlit as st
import pandas as pd

from src.common import *
from src.fileupload import *

params = page_setup()

# Make sure "selected-files" is in session state
if "selected-files" not in st.session_state:
    st.session_state["selected-files"] = params["selected-files"]

st.title("File Upload")

tabs = ["File Upload", "Example Data"]

tabs = st.tabs(tabs)

with tabs[0]:
    dir = {}
    with st.form("Know_RNA", clear_on_submit=True):
        file1 = st.file_uploader("LncRNA_known")
        file2 = st.file_uploader("mRNA_known")
        cols = st.columns(3)
        if cols[1].form_submit_button("Add files to workspace", type="primary"):
            lncRNA_known(file1)
            mRNA_known(file2)  
        try:
           dir['Known_lncRNA'] = file1.name
           dir['Known_mRNA'] = file2.name
        except AttributeError:
           pass
    with st.form("Investigated_sequences", clear_on_submit=True):
        file1 = st.file_uploader("Investigated_sequences")
        cols = st.columns(3)
        if cols[1].form_submit_button("Add files to workspace", type="primary"):
            investigated_sequences(file1)
        try:
           dir['Investigated_sequences'] = file1.name
        except AttributeError:
           pass
# Example mzML files
with tabs[1]:
    st.markdown("Example data set of bacterial cytosolic fractions. Bacillus subtilis cultures were treated with the antibiotic fosfomycin, which inhibits a step in the biosynthesis of petidoglycan (bacterial cell wall). The major accumulation product is UDP-GlcNAc.")
    cols = st.columns(3)
    if cols[1].button("Load Example Data", type="primary"):
        load_example_mzML_files()


if any(Path(upload_files).iterdir()):
    v_space(2)
    # Display all mzML files currently in workspace
    print(st.session_state.workspace)
    df = pd.DataFrame.from_dict(dir, orient = 'index')
    df = df.rename(columns={0: 'Uploaded files'})
    df.to_csv(str(st.session_state.workspace) + '/sample.csv', index=True, sep='\t', index_label='types')
    st.markdown("##### Uploaded files:")
    show_table(df)
    v_space(1)
    # Remove files
    with st.expander("🗑️ Remove files"):
        to_remove = st.multiselect("select files",
                                   options=[f.stem for f in sorted(upload_files.iterdir())])
        c1, c2 = st.columns(2)
        if c2.button("Remove **selected**", type="primary", disabled=not any(to_remove)):
            remove_selected_files(to_remove)
            st.experimental_rerun()

        if c1.button("⚠️ Remove **all**", disabled=not any(upload_files.iterdir())):
            remove_all_files()
            st.experimental_rerun()

save_params(params)
