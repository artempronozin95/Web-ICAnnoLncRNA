import pandas as pd
from Bio import SeqIO
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px

def pie_plot(file1, file2):
    lncfinder = pd.read_csv(file1 + "/lncFinder.csv", sep=',', header=None)
    record_dict = SeqIO.to_dict(SeqIO.parse(file2, "fasta"))
    lncrna = open(str(st.session_state.workspace) + '/Upload_files/LncRNA_prediction/' + 'LncRNA.fasta', 'w', encoding='utf-8')
    coding = open(str(st.session_state.workspace) + '/Upload_files/LncRNA_prediction/' + 'Coding.fasta', 'w', encoding='utf-8')

    lncfinder = lncfinder[[0,1]]
    finder_cod = lncfinder[lncfinder[1].str.contains("NonCoding")==False]
    finder_non = lncfinder[lncfinder[1].str.contains("NonCoding")]
    cod_noncod = [len(finder_cod[1]), len(finder_non[1])]
    for w in finder_non[0]:
        try:
            print(record_dict[w].format('fasta'), end='', file=lncrna)
        except KeyError:
            continue

    for w in finder_cod[0]:
        try:
            print(record_dict[w].format('fasta'), end='', file=coding)
        except KeyError:
            continue


    lables = ['Coding', 'NonCoding']
    fig = px.pie(values=cod_noncod, names=lables, height=300, width=500)
    fig.update_layout(margin=dict(l=50, r=20, t=3, b=0), )
    st.plotly_chart(fig)

