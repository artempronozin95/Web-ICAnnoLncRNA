import pandas as pd
from collections import Counter
from scipy import stats
import statistics
import numpy as np
import plotly.express as px
from pathlib import Path
import streamlit as st
import os
from Bio import SeqIO
import subprocess

from src.common import reset_directory

upload_files: Path = Path(st.session_state.workspace, "Upload_files")

def to_percent(y, position, n):
    s = str(round(100 * y / n, 3))
    return s + '%'

def coord(x):
    Row_list_query = []
    for index, rows in x.iterrows():
        my_list = [rows[1], rows[2]]
        Row_list_query.append(my_list)
    return(Row_list_query)

def freedman_diaconis(data, returnas="width"):
    """
    Use Freedman Diaconis rule to compute optimal histogram bin width.
    ``returnas`` can be one of "width" or "bins", indicating whether
    the bin width or number of bins should be returned respectively.

    Parameters
    ----------
    data: np.ndarray
        One-dimensional array.

    returnas: {"width", "bins"}
        If "width", return the estimated width for each histogram bin.
        If "bins", return the number of bins suggested by rule.
    """
    data = np.asarray(data, dtype=np.float_)
    IQR  = stats.iqr(data, rng=(25, 75), scale="raw", nan_policy="omit")
    N= data.size
    bw   = (2 * IQR) / np.power(N, 1/3)

    if returnas=="width":
        result = bw
    else:
        datmin, datmax = data.min(), data.max()
        datrng=datmax - datmin
        result = int((datrng/bw) + 1)*3
    return result

def gff_prepare(x,y):
    x = x[x['qry_gene_id'].str.contains("path1")]
    name_x = x['qry_gene_id'].str.rsplit('.', expand=True)
    x[['name', 'seq']] = name_x[[0, 1]]
    x = x[x['name'].isin(y[10])]
    x = x.drop_duplicates(subset=['name'])
    x['group'] = np.nan
    x['group'][x['class_code'].isin(['x'])] = "exon antisense"
    x['group'][x['class_code'].isin(['i'])] = "intron"
    x['group'][x['class_code'].isin(['u'])] = "intergenic"
    x.dropna(subset=["group"], inplace=True)
    return x

def plots(query, tmap, small_intron, statistic):
    query = pd.read_csv(query, sep='\t', header=None)
    query = query[pd.to_numeric(query[0], errors='coerce').notnull()]
    tmap = pd.read_csv(tmap, sep='\t')
    if os.stat(small_intron).st_size == 0:
        pass
    else:
        small_intron = pd.read_csv(small_intron, sep='\t', header=None)
        query = query[~query[10].isin(small_intron[0])]
    statistic = open(statistic, 'w', encoding='utf-8')
    tmap = gff_prepare(tmap, query)
    counts = tmap['group'].value_counts()
    counts2 = counts.reset_index()
    fig = px.bar(counts2, x=counts2['index'], y=counts2['group'],
    labels = {'group': 'LncRNA transcripts number', 'index': 'LncRNA classes'}, height = 400)
    st.plotly_chart(fig)


# exon structure of lncRNA
    chart = query[query[7].str.contains("exon")]
    spl = chart[3].str.split('.',2, expand=True)
    size = spl.pivot_table(index = [0], aggfunc ='size').to_dict()
    res = Counter(size.values())
    dict={}
    for k, v in res.items():
        dict[k] = v/len(size)*100
    fig2 = px.bar(x=list(dict.keys()), y=dict.values(),
                 labels={'x': 'Number of exon', 'y': 'Proportion of lncRNAs (%)'}, height=400)
    fig2.update_xaxes(tickvals = np.sort(range(max(list(dict.keys())) + 1)))
    st.plotly_chart(fig2)

#exon size chart
    chart['exon_length'] = (chart[2] - chart[1])
    bin = freedman_diaconis(chart['exon_length'], returnas="bins")
    fig3 = px.histogram(chart, x='exon_length', log_x=True, nbins=bin, histnorm='probability density')
    st.plotly_chart(fig3)

    print('Mean_exon_length','\t', round(statistics.mean(chart['exon_length'])), file=statistic)
    print('Max_exon_length','\t', max(chart['exon_length']), file=statistic)
    print('Min_exon_length','\t', min(chart['exon_length']), file=statistic)
    print('Median_exon_length','\t', statistics.median(chart['exon_length']), file=statistic)

# intron size chart
    chart['name'] = spl[0]
    gr = chart[[1,2,'name']].groupby(by='name')
    intron_length = []
    for key, item in gr:
        if 0 < len(item['name']) - 2:
            n = 0
            intron = abs(item[2].iloc[n] - item[1].iloc[n + 1])
            if intron < 60:
                pass
            elif intron > 100000:
                pass
            else:
                n = n + 1
                intron_length.append(intron)
        else:
            intron_length.append(0)

    intron_length = list(filter(lambda x: x != 0, intron_length))
    if len(intron_length) == 0:
        pass
    else:

        bin_int = freedman_diaconis(intron_length, returnas="bins")
        fig4 = px.histogram(x=intron_length, log_x=True, nbins=bin_int, histnorm='probability density')
        st.plotly_chart(fig4)

        print('Mean_intron_length' ,'\t', statistics.mean(intron_length), file=statistic)
        print('Max_intron_length' ,'\t', max(intron_length), file=statistic)
        print('Min_intron_length' ,'\t', min(intron_length), file=statistic)
        print('Median_intron_length' ,'\t', statistics.median(intron_length), file=statistic)


# choose only query transcripts
    query = query[query[3].str.contains("path1")]
    query['length'] = abs(query[1] - query[2])

    print('Mean transcript length' ,'\t', round(statistics.mean(query['length'])), file=statistic)

# chart lcnRNA distribution across chromosome
    labels, counts = np.unique(query[0].astype(int), return_counts=True)
    fig5 = px.bar(x=labels, y=counts,
    labels = {'x': 'Chromosome of organism', 'y': 'Number of lncRNAs'}, height = 400)
    st.plotly_chart(fig5)

def fasta(x, y, z, r):
    query = x['name']
    for w in query:
        try:
            old_id = y[y['new_id'].isin([w])]
            print(r[str(old_id['old_id'].to_list()[0])].format('fasta'), end='', file=z)
        except KeyError:
            continue

def data_base(x):
    x_n = x.rsplit('/', 1)[1]
    x_n = x_n.rsplit('.', 1)[0]
    index_path = os.path.join('data/reference/data_index/', x_n)
    return index_path

def blast(x,y,evalue, max_target, identity,outfmt):

    path = x.rsplit('/', 1)[0]
    outfmt = os.path.join(path, 'blast' + '.outfmt6')
    aling = 'blastn -query {q} -db {dbw} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"  -evalue {evalue}  -max_target_seqs {max_target} -perc_identity {identity} -out {outfmt}'. format(q=x, dbw=y, evalue=evalue, max_target=max_target, identity=identity, outfmt=outfmt)
    subprocess.call(aling, shell=True)

def aligment(query, tmap, record_dict, old_new, lncrna, evalue, max_target, identity):
    query = pd.read_csv(query, sep='\t', header=None)
    query = query[pd.to_numeric(query[0], errors='coerce').notnull()]
    tmap = pd.read_csv(tmap, sep='\t')
    record_dict = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + record_dict['Uploaded files'][2])
    record_dict = SeqIO.to_dict(SeqIO.parse(record_dict, "fasta"))
    old_new = pd.read_csv(old_new, sep='\t')
    lncrna_file = open(lncrna + 'True_lnc.fasta', 'w', encoding='utf-8')
    tmap = gff_prepare(tmap, query)
    fasta(tmap, old_new, lncrna_file, record_dict)
    true_lnc = os.path.join(str(st.session_state.workspace) + '/Upload_files/new_lncrna/True_lnc.fasta')
    out = os.path.join(lncrna, 'blast.outfmt6')
    database = Path("database/LncAPDB")
    blast(true_lnc, database, evalue[0], max_target[0], identity[0], out)

def add_to_selected_files(filename: str):
    if filename not in st.session_state["ref-files"]:
        st.session_state["ref-files"].append(filename)

def SRX(file):
    if file is not None:
        if file.name not in upload_files.iterdir():
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Add library information")
        return

def value_prepare(x, z, df):
    name = x[0].str.rsplit('_', 1, expand=True)
    x['srx'] = name[1]
    x = x[[0, 'srx']]
    same = pd.merge(df, x, how='inner', left_on=[1], right_on=['srx'])
    pr = same[2].value_counts().reset_index()
    pr = pr.rename(columns={2: z}, inplace=False)
    return (pr, same)

def tissue(blast, old_new, df, org, org_id, out):
    blast = pd.read_csv(blast, sep='\t', header=None)
    LncAPDB = pd.read_csv(Path("database/index_and_newindex.csv"), sep=' ', header=None)
    old_new = pd.read_csv(old_new, sep='\t', header=None, skiprows=1)
    LncAPDB = LncAPDB.drop_duplicates(subset=[1])
    blast_100 = blast[blast[2].isin([100.000])]
    LncAPDB_vs_blast = pd.merge(blast_100, LncAPDB, how='inner', left_on=[1], right_on=[1])
    LncAPDB_vs_blast = LncAPDB_vs_blast[['0_x', 1, '2_x', 10, '0_y', '2_y']]
    LncAPDB_vs_blast.rename(columns={'0_x': 'transcript_id', 1: 'library_id', '2_x': 'percent_identity', 10: 'e_value',
                                     '0_y': 'id_of_database', '2_y': 'database'}, inplace=True)
    out_LncAPDB_vs_blast = os.path.join(out + '/new_lncrna/LncAPDB_vs_blast.csv')
    LncAPDB_vs_blast.to_csv(out_LncAPDB_vs_blast, sep='\t', index=False)
    st.dataframe(LncAPDB_vs_blast)

    # tissue analysis
    lncFinder = pd.read_csv(os.path.join(str(st.session_state.workspace) + '/Upload_files/LncRNA_prediction/lncFinder.csv'), sep=',', header=None)
    coding = lncFinder[~lncFinder[1].isin(['NonCoding'])]

    nonconserved = blast[blast[1].str.contains(org[0])]
    conserved = blast[~blast[1].str.contains(org[0])]

    SRX = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + df)
    df = pd.read_csv(SRX, sep='\t', header=None)
    df = df[df[0].isin([str(org_id)])]
    df_tis = df[2].value_counts().reset_index()
    coding = old_new[old_new[1].isin(coding[0])]
    transc, same = value_prepare(conserved, 'Conserved', df)
    transc_non, same_non = value_prepare(nonconserved, 'Nonconserved',df)
    transc_cod, same_cod = value_prepare(coding, 'Coding', df)

    full_table = pd.merge(transc, transc_cod, how='outer', left_on=['index'], right_on=['index'])
    full_table = pd.merge(full_table, transc_non, how='outer', left_on=['index'], right_on=['index'])
    full_table = pd.merge(full_table, df_tis, how='inner', left_on=['index'], right_on=['index'])
    full_table = full_table.fillna(0)

    full_table['con'] = full_table['Conserved'] / full_table[2]
    full_table['noncon'] = full_table['Nonconserved'] / full_table[2]
    full_table['cod'] = full_table['Coding'] / full_table[2]

    full_table['lncRNA conserved'] = full_table['con'] / len(same[2])
    full_table['lncRNA nonconseved'] = full_table['noncon'] / len(same_non[2])
    full_table['mRNA'] = full_table['cod'] / len(same_cod[2])

    full_table = full_table[['index', 'lncRNA conserved', 'lncRNA nonconseved', 'mRNA']]

    full_table = full_table.set_index('index')
    fig = px.imshow(full_table, aspect="auto")
    st.plotly_chart(fig)

