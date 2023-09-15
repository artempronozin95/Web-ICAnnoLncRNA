import streamlit as st
from streamlit_plotly_events import plotly_events
import subprocess
import shutil
from pathlib import Path

from src.common import *
from src.test_train import *
from src.TP_FN import *
from src.rename import *
from src.take_noncoding import *

params = page_setup()
st.title("LncRNA prediction")

upload_files: Path = Path(st.session_state.workspace, "Upload_files")

selected_file = pd.read_csv(str(st.session_state.workspace) + '/sample.csv', sep='\t')
st.markdown("##### Uploaded files:")
st.dataframe(selected_file, hide_index=True, use_container_width=True)

with st.form("LncRNA_prediction"):
	c1, c2 = st.columns(2)
	c1.radio(
            "Structure type",
            ["DNA", "SS"],
            ["DNA", "SS"].index(params["structure_type"]),
            key="structure_type",
            help="Need to choose whether to use secondary structure in model building or not. Choose between DNA or SS (secondary structure).",)
	
	c1, c2, c3 = st.columns(3)
	if c2.form_submit_button("Run prediction", type="primary"):
	        path = os.path.join(str(st.session_state.workspace) + '/Upload_files/test_train') 
	        
	        if os.path.exists(path) and os.path.isdir(path):
	            shutil.rmtree(path)
	        else:
	            pass  
	        
	        known_lncRNA = selected_file[selected_file['types'].str.contains("Known_lncRNA")]
	        known_mRNA = selected_file[selected_file['types'].str.contains("Known_mRNA")]
	        test_train(known_lncRNA['Uploaded files'][1], known_mRNA['Uploaded files'][0])  
	        structure = params['structure_type']
	        with st.spinner("Models building..."):
	            process = subprocess.Popen(["Rscript", "./src/model_lncFind.r", path, structure], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	            result = process.communicate()
	        st.success('Models building done', icon="✅")
		
	        with st.spinner("Best model is..."):
	            df, best = best_model(path)
	            st.markdown("##### Best model accuracy:")
	            st.dataframe(df, hide_index=True, use_container_width=True)
	            model = os.path.join(best + 'model.rds')

	        results = os.path.join(str(st.session_state.workspace) + '/Upload_files/LncRNA_prediction')
	        os.makedirs(results, exist_ok=True)
	        Investigated_sequences = selected_file[selected_file['types'].str.contains("Investigated_sequences")]
	        investigated_file = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + str(Investigated_sequences['Uploaded files'][2]))
	        corrected_file = os.path.join(str(st.session_state.workspace) + '/Upload_files/filtered_seq.fasta')
	        new_vs_old_id = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_vs_old_id.tsv')
	        with st.spinner("Investigated file preparing..."):
	            prepare_seq(investigated_file, corrected_file, new_vs_old_id)	
			
	        with st.spinner("LncRNA prediction..."):
	            process2 = subprocess.Popen(["Rscript", "./src/lncFind.r", best, structure, corrected_file, results])
	            result2 = process2.communicate()
	            pie_plot(results, corrected_file)
	        st.success('LncRNA prediction complete', icon="✅")



try:
   st.download_button(label="Download Model", data=model, file_name='model.rds', mime="model/rds")
except NameError:
   pass
