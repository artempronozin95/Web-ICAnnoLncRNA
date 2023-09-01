import numpy as np
from sklearn.model_selection import train_test_split
from Bio import SeqIO
import os
import streamlit as st
from pathlib import Path
from src.common import reset_directory

def test_train(file1, file2):
    with st.spinner("Files parsing..."):
        lnc = SeqIO.to_dict(SeqIO.parse(str(st.session_state.workspace) + '/Upload_files/' + str(file1), "fasta"))
        cds = SeqIO.to_dict(SeqIO.parse(str(st.session_state.workspace) + '/Upload_files/' + str(file2), "fasta"))
    with st.spinner("Building training/test data..."):
        directory = os.path.join(str(st.session_state.workspace) + '/Upload_files/test_train/')
        os.makedirs(directory, exist_ok=True)
        train_lnc = open(directory + 'train_lnc.fasta', 'w', encoding='utf-8')
        test_lnc = len(lnc.keys())*0.25
        test_mrna = test_lnc*2
        n=0
        keys = list(cds.keys())
        pers_lnc=test_lnc/len((list(lnc.keys())))
        Y_train, Y_test = train_test_split(list(lnc.keys()), test_size=pers_lnc, shuffle=True)
        if len(Y_train)*2 <= len(cds.keys())/5:
            train_mrna = len(Y_train)*2
        else:
            train_mrna = len(cds.keys())/5.2
        for w in Y_train:
            print(lnc[w].format('fasta'), end='', file=train_lnc)
        while n <= 4:
            os.mkdir(directory + str(n))
            train_cds = open(directory + str(n) + '/train_mrna.fasta', 'w', encoding='utf-8')
            test_file = open(directory + str(n) + '/test.fasta', 'w', encoding='utf-8')
            compare = open(directory + str(n) + '/compare.csv', 'w', encoding='utf-8')
            pers_mrna=train_mrna/(len(keys))
            X_train, X_test = train_test_split(keys, test_size=pers_mrna, shuffle=True)
            print('train_set',len(X_train))
            print('test_set',len(X_test))
            if len(X_test) > test_mrna:
    	        pers_mrna_test = test_mrna/(len(X_test))
            else:
                true_per = len(X_test)*0.25
                pers_mrna_test = true_per/len(X_test)
            pers_mrna_test = test_mrna/(len(X_test))
            train, test = train_test_split(X_test, test_size=pers_mrna_test, shuffle=True)
            for w in test:
                print(cds[w].format('fasta'), end='', file=test_file)
                print(cds[w].id , 'mrna', sep='\t', file=compare)
            for w in train:
                print(cds[w].format('fasta'), end='', file=train_cds)
            for w in Y_test:
                print(lnc[w].format('fasta'), end='', file=test_file)
                print(lnc[w].id, 'lnc', sep='\t', file=compare)
            keys = [w for w in keys if w not in X_test]
            print('new_set',len(keys))
            n=n+1
