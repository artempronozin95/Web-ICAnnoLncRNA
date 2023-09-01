from pathlib import Path
import pandas as pd

import streamlit as st
from streamlit_plotly_events import plotly_events
import subprocess

from src.common import *
from src.statistic import *

params = page_setup()
st.title("Statistic for predicted lncRNA")

upload_files: Path = Path(st.session_state.workspace, "Upload_files")
selected_file = pd.read_csv(str(st.session_state.workspace) + '/sample.csv', sep='\t')

query = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/True_lnc.bed')
tmap = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/gffcmp_first.filter_map.bed.tmap')
introns = os.path.join(str(st.session_state.workspace) + '/Upload_files/GMAP/itron_small.txt')
out = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/stats.txt')

plots(query, tmap, introns, out)

df = pd.read_csv(out, sep='\t', header=None)
df1 = df.rename(columns={0: "Structural Characteristic", 1: "Value"})
st.dataframe(df1, hide_index=True, use_container_width=True)

record_dict = selected_file[selected_file['types'].str.contains("Investigated_sequences")]
old_new = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_vs_old_id.tsv')
lncrna = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/')

with st.form("Homologs search", clear_on_submit=True):
         c1, c2, c3 = st.columns(3)
         evalue = c1.text_input("evalue", "")
         max_target = c2.text_input("max_target_seq", "")
         identity = c3.text_input("perc_identity", "")
         if c2.form_submit_button("Run alignment", type="primary"):
            aligment(query, tmap, record_dict, old_new, lncrna, evalue, max_target, identity)
         st.success('Alignment done', icon="✅")
         
with st.form("Library information", clear_on_submit=True):
     file = st.file_uploader("Library information")
     cols = st.columns(3)
     dir = {}
     if cols[1].form_submit_button("Add files to workspace", type="primary"):
        SRX(file)
     try:
        dir['Library information'] = file.name             
     except AttributeError:
        pass
    
     df2 = pd.DataFrame.from_dict(dir, orient = 'index')
     df2 = df2.rename(columns={0: 'Uploaded files'})
     st.markdown("##### Uploaded files:")
     show_table(df2)

with st.form("Tissue specificity", clear_on_submit=True):
     c1, c2 = st.columns(2)
     org = c1.text_input("Organism", "")
     org_id = c2.text_input("Organism_ID", "")
     blast = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/blast.outfmt6')
     out = os.path.join(str(st.session_state.workspace) + '/Upload_files')
     c1, c2, c3 = st.columns(3)
     if c2.form_submit_button("Tissue specificity analysis", type="primary"):
        tissue(blast, old_new, df2['Uploaded files'][0], org, org_id, out)
     
     
