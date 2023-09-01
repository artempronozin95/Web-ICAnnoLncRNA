from pathlib import Path
import streamlit as st
from src.common import reset_directory

upload_files: Path = Path(st.session_state.workspace, "Upload_files")

def add_to_selected_files(filename: str):

    if filename not in st.session_state["selected-files"]:
        st.session_state["selected-files"].append(filename)


def lncRNA_known(file):
    if file is not None:
        if file.name not in upload_files.iterdir() and file.name.endswith("fasta"):
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Add first file")
        return

def mRNA_known(file):
    if file is not None:
        if file.name not in upload_files.iterdir() and file.name.endswith("fasta"):
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Add second file")
        return

def investigated_sequences(file):
    if file is not None:
        if file.name not in upload_files.iterdir() and file.name.endswith("fasta"):
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Upload a file first")
        return


def remove_selected_files(to_remove: list[str]) -> None:
    # remove all given files from mzML workspace directory and selected files
    for f in to_remove:
        Path(upload_files, f+".fasta").unlink()
        st.session_state["selected-files"].remove(f)
    st.success("Selected files removed!")


def remove_all_files() -> None:
    # reset (delete and re-create) mzML directory in workspace
    reset_directory(upload_files)
    # reset selected mzML list
    st.session_state["selected-files"] = []
    st.success("All files removed!")
     




