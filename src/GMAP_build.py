import os
import subprocess
import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
from pybedtools import BedTool

from src.common import reset_directory

upload_files: Path = Path(st.session_state.workspace, "Upload_files")

def add_to_selected_files(filename: str):
    if filename not in st.session_state["ref-files"]:
        st.session_state["ref-files"].append(filename)

def ref_genome(file):
    if file is not None:
        if file.name not in upload_files.iterdir():
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Add genome file")
        return

def annotation(file):
    if file is not None:
        if file.name not in upload_files.iterdir():
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Add annotation file")
        return

def TE(file):
    if file is not None:
        if file.name not in upload_files.iterdir():
            with open(Path(upload_files, file.name), "wb") as fh:
                fh.write(file.getbuffer())
        add_to_selected_files(Path(file.name).stem)
    if file is None:
        st.warning("Add transposable elements coordinates file")
        return

def build(reference):
    ref = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + reference)
    size = round(os.path.getsize(ref) / (1024 * 1024))
    if size <= 64:
        k = 12
    elif 64 < size <= 256:
        k = 13
    elif 256 < size <= 1000:
        k = 14
    elif 1000 < size <= 4000:
        k = 15

    gmap_build = 'gmap_build'
    path, base = os.path.split(ref)
    ref_label = os.path.splitext(base)[0]
    command = '{gmap_build} -D {tmp_dir} -d {ref_index_name} -k {kmer_value} {reference}'.\
            format(gmap_build=gmap_build, tmp_dir=path, ref_index_name=ref_label, reference=ref, kmer_value=k)
    subprocess.call(command, shell=True)

def aling(reference, transcripts, min, max, lg, out, gff):
    ref = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + reference)
    gff = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + gff)
    if int(lg) > 2**32:
       gmap_run = 'gmapl'
       print('large genome use gmapl')
    else:
       gmap_run = 'gmap'
    path, base = os.path.split(ref)
    ref_label = os.path.splitext(base)[0]
    out_res = os.path.join(out, 'gmap_mapping.gff')
    gmap_run_logger_err_path = os.path.join(out, gmap_run + '.out.log')
    command = '{gmap} -D {tmp_dir} -d {ref_index_name} {transcripts} --min-intronlength={min_intron_length} --intronlength={intron_length} --cross-species ' \
                  '--format=gff3_gene --split-large-introns --npaths=1 > {alignment_out} 2>> {log_out_2}'.\
            format(gmap=gmap_run, tmp_dir=path, ref_index_name=ref_label, transcripts=transcripts,
                   alignment_out=out_res, log_out_2=gmap_run_logger_err_path, min_intron_length=[0], intron_length=max[0])
    bed = out_res.rsplit('.', 1)
    bed = bed[0]
    ref_bed = gff.rsplit('.', 1)
    ref_bed = ref_bed[0]
    command1 = '{gff2bed} < {query} > {query_out}'. format(gff2bed = 'gff2bed', query = out_res, query_out = bed + '.bed')
    command2 = '{gff2bed} < {reference} > {reference_out}'. format(gff2bed = 'gff2bed', reference = gff, reference_out = ref_bed + '.bed')
    subprocess.call(command, shell=True)
    subprocess.call(command1, shell=True)
    subprocess.call(command2, shell=True)


def filter_introns(gmap, out):
    out_res = os.path.join(out, 'filter_map.bed')
    gmap_align = pd.read_csv(gmap, sep='\t', header=None)
    name = gmap_align[9].str.split(';', 4, expand=True)
    name = name[1].str.split('=', 2, expand=True)
    gmap_align['name'] = name[1]
    exon = gmap_align[gmap_align[3].str.contains("exon")]
    gr = exon[[0, 1, 2, 'name']].groupby(by='name')
    intron_large = []
    intron_small = []
    intron_data = []
    for key, item in gr:
        n = 0
        if 0 < len(item['name']) - 2:
            intron = abs(item[2].iloc[n] - item[1].iloc[n + 1])
            intron_data.append(
                [item[0].iloc[n], item[2].iloc[n], item[1].iloc[n + 1], key + '_intron' + str(n + 1), intron])
            n = n + 1
            if intron < 60:
                intron_small.append(key)
            if intron > 20000:
                intron_large.append(key)
            else:
                pass
        else:
            pass

    its = pd.DataFrame(intron_small)
    its.to_csv(os.path.join(out, 'itron_small.txt'), index=False, header=None)
    gmap_align = gmap_align[~gmap_align['name'].isin(intron_large)]
    gmap_align['length'] = gmap_align[2] - gmap_align[1]
    gmap_align = gmap_align[gmap_align['length'] < 5000]
    gmap_align_big = gmap_align[gmap_align['length'] > 5000]
    gmap_align = gmap_align[~gmap_align['name'].isin(gmap_align_big['name'])]
    gmap_align = gmap_align[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    gmap_align.to_csv(out_res, index=False, sep='\t', header=None)

def gff_compare(ref, map, out, index):
    ref = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + ref)
    command3 = '{gffcompare} -r {mapping} {query} -o {out}'.format(gffcompare='gffcompare', mapping=ref, query=map, out=os.path.join(out, 'gffcmp_' + index))
    subprocess.call(command3, shell=True)

def clean(x):
    split = x.str.rsplit('.', 2, expand=True)
    split = split[0]
    return (split)

def prepare_table(tmap,gmap):
    query = clean(tmap['qry_gene_id'])
    gmap['name'] = clean(gmap[3])
    query = query.drop_duplicates()
    # take lncRNA from alignment file
    align = gmap[gmap['name'].isin(query)]
    align = align[align[7].str.contains("gene")]
    return align, query

def loci(tmap, gmap, out):
    tmap = pd.read_csv(tmap, sep='\t')
    gmap_align = pd.read_csv(gmap, sep='\t', header=None)
# take only lncRNA classes from tmap file
    tmap['group'] = np.nan
    tmap['group'][tmap['class_code'].isin(['x'])] = "exon antisense"
    tmap['group'][tmap['class_code'].isin(['i'])] = "intron"
    tmap['group'][tmap['class_code'].isin(['u'])] = "intergenic"
    tmap.dropna(subset = ["group"], inplace=True)
    tmap = tmap[tmap['qry_gene_id'].str.contains("path1")]
# prepare table
    before_loci, q = prepare_table(tmap, gmap_align)
    before_loci.to_csv(os.path.join(out + 'before_loci.bed'), sep='\t', index=False, header=None)
# lncRNA into loci
    genes = BedTool(os.path.join(out + 'before_loci.bed'))
    loci = genes.merge(s=True, c=[4,5,6], o=['collapse','mean','distinct'])
    loci_df = pd.read_table(loci.fn, sep='\t', names=[0,1,2,3,4,5])
    loci_df = loci_df[(loci_df[2] - loci_df[1])>200]
    loci_df['numb'] = list(range(0, len(loci_df[3])))
    loci_df['loc'] = 'LOC'
    loci_df['loci'] = loci_df['loc'] + '_' + loci_df['numb'].astype(str)
    loci_df = loci_df[[0,1,2,'loci',4,5,3]]
    loci_df.to_csv(os.path.join(out + 'after_loci.bed'), sep='\t', index=False, header=None)
    out_dir = os.path.join(str(st.session_state.workspace) + '/Upload_files/gff_second/')
    os.makedirs(out_dir, exist_ok=True)
    gff_compare(os.path.join('loci/' + 'after_loci.bed'), gmap, out_dir, 'second')

def TE_filter(tmap, gmap, TE, out, true_lnc):
    TE = os.path.join(str(st.session_state.workspace) + '/Upload_files/' + TE)
    tmap = pd.read_csv(tmap, sep='\t')
    gmap_align_prime = pd.read_csv(gmap, sep='\t', header=None)
# take only lncRNA classes from tmap file
    tmap['group'] = np.nan
    tmap['group'][tmap['class_code'].isin(['='])] = 'true lncrna'
    tmap.dropna(subset=["group"], inplace=True)
    tmap = tmap[tmap['qry_gene_id'].str.contains("path1")]
# prepare table
    prepare_tab, q = prepare_table(tmap, gmap_align_prime)
    prepare_tab.to_csv(os.path.join(out + 'lnc_after_loci.bed'), sep='\t', index=False, header=None)
# BedTool intersect
    genes = BedTool(os.path.join(out + 'lnc_after_loci.bed'))
    loci = genes.intersect(wa=True, wb=True, b=TE)
    loci_df = pd.read_table(loci.fn, sep='\t',
                                names=['chr', 'lncRNA_start', 'lncRNA_end', 'lncRNA', 4, 'lncRNA_strand', 6, 7, 8, 9,
                                       10, 'chr_TE', 'TE_start', 'TE_end'])
    loci_df = loci_df[
            ['chr', 'lncRNA_start', 'lncRNA_end', 'lncRNA', 'lncRNA_strand', 'chr_TE', 'TE_start', 'TE_end']]
    loci_df = loci_df.drop_duplicates(subset=['lncRNA'])
    loci_name = loci_df['lncRNA'].str.split('.', 4, expand=True)
    true_lncrna = gmap_align_prime[~gmap_align_prime['name'].isin(loci_name[0])]
    true_lncrna.to_csv(os.path.join(true_lnc + 'True_lnc.bed'), sep='\t', index=False, header=None)

def true_lnc(tmap, gmap, true_lnc):
    tmap = pd.read_csv(tmap, sep='\t')
    gmap_align_prime = pd.read_csv(gmap, sep='\t', header=None)
    # take only lncRNA classes from tmap file
    tmap['group'] = np.nan
    tmap['group'][tmap['class_code'].isin(['='])] = 'true lncrna'
    tmap.dropna(subset=["group"], inplace=True)
    tmap = tmap[tmap['qry_gene_id'].str.contains("path1")]
    # prepare table
    prepare_tab, q = prepare_table(tmap, gmap_align_prime)
    gmap_align_prime = gmap_align_prime[gmap_align_prime['name'].isin(q)]
    gmap_align_prime.to_csv(os.path.join(true_lnc + 'True_lnc.bed'), sep='\t', index=False, header=None)


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


