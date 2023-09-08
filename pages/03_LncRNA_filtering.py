from pathlib import Path
from os import listdir
from pathlib import PurePath

import streamlit as st
from streamlit_plotly_events import plotly_events
import subprocess

from src.common import *
from src.GMAP_build import *

params = page_setup()
st.title("LncRNA filtering")

upload_files: Path = Path(st.session_state.workspace, "Upload_files")
os.makedirs(os.path.join(str(st.session_state.workspace) + '/Upload_files/genome/'), exist_ok=True)
genome_id: Path = Path(st.session_state.workspace, "Upload_files/genome")

if "ref-files" not in st.session_state:
    st.session_state["ref-files"] = params["ref-files"]


with st.expander("Settings", expanded=True):
    c1, c2 = st.columns(2)
    c1.radio(
            "Genome index",
            ["Genome index built", "No genome index"],
            ["Genome index built", "No genome index"].index(params["index"]),
            key="index",
            help="",)

    dir = {}
    if params['index'] == "No genome index":
    	with st.form("Reference_genome", clear_on_submit=True):
            file1 = st.file_uploader("Reference genome")
            file2 = st.file_uploader("Reference genome annotation (gff, gtf, gff3)")
            cols = st.columns(3)
            if cols[1].form_submit_button("Add files to workspace", type="primary"):
               ref_genome(file1)
               annotation(file2)
            try:
               dir['Reference genome'] = file1.name
               dir['Reference genome annotation'] = file2.name                 
            except AttributeError:
               pass
            df = pd.DataFrame.from_dict(dir, orient = 'index')
            df = df.rename(columns={0: 'Uploaded files'})
            df.to_csv(str(st.session_state.workspace) + '/sample.csv', index=True, sep='\t', index_label='types')
            st.markdown("##### Uploaded files:")
            show_table(df)
    	with st.form("Index_genome", clear_on_submit=True):
    	    cols = st.columns(3)
    	    if cols[1].form_submit_button("Run index build", type="primary"):
    	        with st.spinner("Indexing..."):
                     build(df['Uploaded files'][0]) 
    	        st.success('Indexing done', icon="✅")
    
    elif params['index'] == "Genome index built":
       list = listdir(Path("genome"))
       selected_file = st.selectbox("choose organism", [f for f in list],)
       for w in Path("genome", selected_file).glob("*.fa"):
           path_split = PurePath(w).parts
           shutil.copy(w, upload_files)
           add_to_selected_files(w.stem)
           dir['Reference genome'] = path_split[2]
       for w in Path("genome", selected_file).glob("*.gff"):
           path_split = PurePath(w).parts
           shutil.copy(w, upload_files)
           add_to_selected_files(w.stem)
           dir['Reference genome annotation'] = path_split[2] 
       for w in Path("genome", selected_file, "genome").glob("*"):
           path_split = PurePath(w).parts
           shutil.copy(w, genome_id)
           add_to_selected_files(w.stem)
       df = pd.DataFrame.from_dict(dir, orient = 'index')
       df = df.rename(columns={0: 'Uploaded files'})
       df.to_csv(str(st.session_state.workspace) + '/sample.csv', index=True, sep='\t', index_label='types')
       st.markdown("##### Uploaded files:")
       show_table(df)
    
    col1, col2 = st.columns(2)
    col1.radio(
            "Transposable element",
            ["is available", "no available"],
            ["is available", "no available"].index(params["Transposable_element"]),
            key="Transposable_element",
            help="",)
    TE_av = []
    if params['Transposable_element'] == "is available":
       TE_av.append('TE')
       with st.form("Transposable_element", clear_on_submit=True):
              file3 = st.file_uploader("Transposable elements (bed)")
              cols = st.columns(3)
              if cols[1].form_submit_button("Add files to workspace", type="primary"):
                 TE(file3)
              try:
                 dir['Transposable element'] = file3.name             
              except AttributeError:
                 pass
              df1 = pd.DataFrame.from_dict(dir, orient = 'index')
              df1 = df1.rename(columns={0: 'Uploaded files'})
              st.markdown("##### Uploaded files:")
              show_table(df1)
    if params['Transposable_element'] == "no available":
       TE_av.append('no_TE')
    
    with st.form("Mapping", clear_on_submit=True):
         c1, c2, c3 = st.columns(3)
         corrected_file = os.path.join(str(st.session_state.workspace) + '/Upload_files/filtered_seq.fasta')
         min_intron_length = c1.text_input("Min intron length", "")
         max_intron_length = c2.text_input("Max intron length", "")
         genome_length = c3.text_input("Genome length", "")
         results = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/')
         os.makedirs(results, exist_ok=True)
         if c2.form_submit_button("Run mapping", type="primary"):
    	     with st.spinner("Mapping..."):
                  aling(df['Uploaded files'][0], corrected_file, min_intron_length, max_intron_length, genome_length, results, df['Uploaded files'][1]) 
    	     st.success('Mapping done', icon="✅")
    	     mapping = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/gmap_mapping.bed')
    	     with st.spinner("Filter sequence with short introns..."):
                  filter_introns(mapping, results) 
    	     st.success('Filtering done', icon="✅")
    	     filt_mapping = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/filter_map.bed')
    	     out = os.path.join(str(st.session_state.workspace) + '/Upload_files/gff_first/')
    	     os.makedirs(out, exist_ok=True)
    	     with st.spinner("Sequence classification..."):
                  gff_compare(df['Uploaded files'][1], filt_mapping, out, 'first') 
    	     st.success('Sequence classification done', icon="✅")
    	     if TE_av[0] == 'TE':
                tmap = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/gffcmp_first.filter_map.bed.tmap')
                filt_mapping = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/filter_map.bed')
                out_filt = os.path.join(str(st.session_state.workspace) + '/Upload_files/loci/')
                os.makedirs(out_filt, exist_ok=True)
            # delete overrepresented sequences
                loci(tmap, filt_mapping, out_filt)
            # TE filter
                tmap2 = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/gffcmp_second.filter_map.bed.tmap')
                TE_out = os.path.join(str(st.session_state.workspace) + '/Upload_files/loci/')
                True_lnc = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/')
                os.makedirs(True_lnc, exist_ok=True)
                os.makedirs(TE_out, exist_ok=True)
                TE_filter(tmap2, filt_mapping, df1['Uploaded files'][2], TE_out, True_lnc)
    	     if TE_av[0] == 'no_TE':
                tmap = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/gffcmp_first.filter_map.bed.tmap')
                filt_mapping = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/filter_map.bed')
                out_filt = os.path.join(str(st.session_state.workspace) + '/Upload_files/loci/')
                os.makedirs(out_filt, exist_ok=True)
            # delete overrepresented sequences
                loci(tmap, filt_mapping, out_filt)
            # TE filter
                tmap2 = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/gffcmp_second.filter_map.bed.tmap')
                TE_out = os.path.join(str(st.session_state.workspace) + '/Upload_files/loci/')
                True_lnc = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/')
                os.makedirs(True_lnc, exist_ok=True)
                os.makedirs(TE_out, exist_ok=True)
                true_lnc(tmap2, filt_mapping, True_lnc)


if any(Path(upload_files).iterdir()):
    v_space(2)
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
