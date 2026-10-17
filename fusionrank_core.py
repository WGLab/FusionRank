#!/usr/bin/env python3

import json
import os
from html import escape
import pickle
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pybedtools

# ╔══════════════════════════════════════════════════════════════╗
# ║  [0]  PREPROCESSING                                          ║
# ╚══════════════════════════════════════════════════════════════╝

BKP_INPUT_COLS = [
    'gene1_chr', 'gene1_bkp', 'gene1_strand',
    'gene2_chr', 'gene2_bkp', 'gene2_strand',
]

def _has_breakpoint_cols(df):
    """Returns True if all required breakpoint columns are present in df."""
    return all(c in df.columns for c in BKP_INPUT_COLS)

def preprocess_fusion_df(input_file):
    print(f'Preprocessing: {input_file}')
    df = pd.read_csv(input_file, sep='\t')

    has_bkp = _has_breakpoint_cols(df)

    if has_bkp:
        # process remaining fusions
        for col in ['gene1_chr', 'gene2_chr']:
            df[col] = (df[col].astype(str).str.strip()
                       .apply(lambda x: x if x.startswith('chr') else f'chr{x}'))
        for col in ['gene1_bkp', 'gene2_bkp']:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
        for col in ['gene1_strand', 'gene2_strand']:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: x if x in ['+', '-'] else None)
    else:
        print('    ⚠️  Breakpoint columns not found — skipping breakpoint preprocessing and cleaning')

    if 'fusion_name' not in df.columns:
        df['fusion_name'] = df['gene1'] + '::' + df['gene2']

    if has_bkp:
        critical = ['gene1_chr', 'gene1_bkp', 'gene1_strand', 'gene2_chr', 'gene2_bkp', 'gene2_strand']
        n_before = len(df)
        df = df.dropna(subset=critical)
        print(f'    {n_before:,} → {len(df):,} rows after cleaning')
    else:
        print(f'    {len(df):,} rows loaded (no breakpoint cleaning applied)')

    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [1]  CIViC                                                  ║
# ║       → in_civic_exact, in_civic_wildcard                    ║
# ║       → civic_best_evidence_score, civic_has_therapeutic     ║
# ║       → civic_match_score                                    ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_civic(df, refs):
    """
    Adds: in_civic_exact, in_civic_wildcard, civic_best_evidence_score,
          civic_has_therapeutic, civic_match_score
    Used by: relevance score (tier 4 therapeutic check)
    """
    print('CIViC...')
    df = df.copy()
    civic = pd.read_csv(refs['civic_file'], sep='\t')

    civic_exact = civic[civic['is_wildcard_fusion'] == False].copy()
    civic_wc    = civic[civic['is_wildcard_fusion'] == True].copy()

    # fusion pair sets (both orientations)
    civic_pairs = (
        set((civic_exact['gene_a'] + '::' + civic_exact['gene_b']).str.upper()) |
        set((civic_exact['gene_b'] + '::' + civic_exact['gene_a']).str.upper())
    )
    civic_wc_genes = set(
        civic_wc[civic_wc['gene_a'] == 'any_gene']['gene_b'].str.upper().tolist() +
        civic_wc[civic_wc['gene_b'] == 'any_gene']['gene_a'].str.upper().tolist()
    )

    # best evidence per fusion
    emap = {'A':5,'B':4,'C':3,'D':2,'E':1}
    civic_exact['ev_score'] = civic_exact['evidence_level'].map(emap)
    best_ev = (civic_exact.groupby(['gene_a','gene_b'])['ev_score'].max()
               .reset_index().rename(columns={'ev_score':'civic_best_evidence_score'}))
    best_ev['fn'] = best_ev['gene_a'] + '::' + best_ev['gene_b']
    ev_map = dict(zip(best_ev['fn'].str.upper(), best_ev['civic_best_evidence_score']))
    ev_map.update(dict(zip((best_ev['gene_b']+'::'+best_ev['gene_a']).str.upper(),
                           best_ev['civic_best_evidence_score'])))

    # therapeutic
    civic_exact_ther = set(
        (civic_exact[civic_exact['evidence_type']=='Predictive']['gene_a'] + '::' +
         civic_exact[civic_exact['evidence_type']=='Predictive']['gene_b']).str.upper()
    )
    civic_wc_ther_genes = set(
        civic_wc[(civic_wc['is_wildcard_fusion']==True) & (civic_wc['evidence_type']=='Predictive') &
                 (civic_wc['gene_a']=='any_gene')]['gene_b'].str.upper().tolist() +
        civic_wc[(civic_wc['is_wildcard_fusion']==True) & (civic_wc['evidence_type']=='Predictive') &
                 (civic_wc['gene_b']=='any_gene')]['gene_a'].str.upper().tolist()
    )

    # cancer-type dicts
    civic_exact_by_cancer = {}
    for ct in civic_exact['tcga_cancer_type'].unique():
        sub = civic_exact[civic_exact['tcga_cancer_type']==ct]
        civic_exact_by_cancer[ct] = (
            set((sub['gene_a']+'::'+sub['gene_b']).str.upper()) |
            set((sub['gene_b']+'::'+sub['gene_a']).str.upper())
        )
    civic_wc_by_cancer = {}
    for ct in civic_wc['tcga_cancer_type'].unique():
        sub = civic_wc[civic_wc['tcga_cancer_type']==ct]
        civic_wc_by_cancer[ct] = set(
            sub[sub['gene_a']=='any_gene']['gene_b'].str.upper().tolist() +
            sub[sub['gene_b']=='any_gene']['gene_a'].str.upper().tolist()
        )

    fn_up = df['fusion_name'].str.upper()
    g1_up = df['gene1'].str.upper()
    g2_up = df['gene2'].str.upper()

    df['in_civic_exact']           = fn_up.isin(civic_pairs).astype(int)
    df['in_civic_wildcard']        = (g1_up.isin(civic_wc_genes) | g2_up.isin(civic_wc_genes)).astype(int)
    df['civic_best_evidence_score']= fn_up.map(ev_map).fillna(0).astype(int)
    df['civic_has_therapeutic']    = (
        fn_up.isin(civic_exact_ther) |
        g1_up.isin(civic_wc_ther_genes) | g2_up.isin(civic_wc_ther_genes)
    ).astype(int)

    def _civic_match(fn, g1, g2, ct):
        fu = fn.upper(); g1u = g1.upper(); g2u = g2.upper()
        if ct in civic_exact_by_cancer and fu in civic_exact_by_cancer[ct]: return 4
        if ct in civic_wc_by_cancer and (g1u in civic_wc_by_cancer[ct] or g2u in civic_wc_by_cancer[ct]): return 3
        if fu in civic_pairs: return 2
        if g1u in civic_wc_genes or g2u in civic_wc_genes: return 1
        return 0

    df['civic_match_score'] = df.apply(
        lambda r: _civic_match(r['fusion_name'], r['gene1'], r['gene2'], r['cancer_type']), axis=1)

    print(f'    in_civic_exact: {df["in_civic_exact"].sum()} | '
          f'has_therapeutic: {df["civic_has_therapeutic"].sum()}')
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [2]  COSMIC                                                 ║
# ║       → in_cosmic_exact, in_cosmic_same_cancer               ║
# ║       → cosmic_gene_involved, cosmic_match_score             ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_cosmic(df, refs):
    """
    Adds: in_cosmic_exact, in_cosmic_same_cancer,
          cosmic_gene_involved, cosmic_match_score
    Used by: relevance score (tier 3 & tier 2 pathway 1)
    """
    print('COSMIC...')
    df = df.copy()
    cosmic = pd.read_csv(refs['cosmic_file'], sep='\t')

    split = cosmic['fusion_name'].str.split('::')
    cosmic_pairs = (
        set(cosmic['fusion_name'].str.upper()) |
        set((split.str[1] + '::' + split.str[0]).str.upper())
    )
    cosmic_genes = set(split.str[0].str.upper()) | set(split.str[1].str.upper())

    cosmic_by_cancer = {}
    for ct in cosmic['tcga_cancer_type'].unique():
        sub = cosmic[cosmic['tcga_cancer_type']==ct]
        sp  = sub['fusion_name'].str.split('::')
        cosmic_by_cancer[ct] = (
            set(sub['fusion_name'].str.upper()) |
            set((sp.str[1]+'::'+sp.str[0]).str.upper())
        )

    fn_up = df['fusion_name'].str.upper()
    df['in_cosmic_exact']      = fn_up.isin(cosmic_pairs).astype(int)
    df['in_cosmic_same_cancer']= df.apply(
        lambda r: int(r['cancer_type'] in cosmic_by_cancer and
                      r['fusion_name'].upper() in cosmic_by_cancer[r['cancer_type']]), axis=1)
    df['cosmic_gene_involved'] = (
        df['gene1'].str.upper().isin(cosmic_genes) |
        df['gene2'].str.upper().isin(cosmic_genes)
    ).astype(int)
    df['cosmic_match_score']   = df.apply(
        lambda r: (2 if r['cancer_type'] in cosmic_by_cancer and
                        r['fusion_name'].upper() in cosmic_by_cancer[r['cancer_type']]
                   else 1 if r['fusion_name'].upper() in cosmic_pairs else 0), axis=1)

    print(f'    in_cosmic_exact: {df["in_cosmic_exact"].sum()} | '
          f'cosmic_gene_involved: {df["cosmic_gene_involved"].sum()}')
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [3]  MITELMAN                                               ║
# ║       → in_mitelman_exact, in_mitelman_same_cancer           ║
# ║       → mitelman_literature_count, mitelman_same_cancer_count║
# ║       → mitelman_match_score                                 ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_mitelman(df, refs):
    """
    Adds: in_mitelman_exact, in_mitelman_same_cancer,
          mitelman_literature_count, mitelman_same_cancer_count,
          mitelman_match_score
    Used by: relevance score (tier 2 pathway 3)
    """
    print('Mitelman...')
    df = df.copy()
    mit = pd.read_csv(refs['mitelman_file'], sep='\t').rename(columns={'Fusion Gene':'fusion_name'})

    split = mit['fusion_name'].str.split('::')
    mit_pairs = (
        set(mit['fusion_name'].str.upper()) |
        set((split.str[1]+'::'+split.str[0]).str.upper())
    )
    max_count = (mit.groupby('fusion_name')['Count'].max().reset_index()
                 .rename(columns={'Count':'mitelman_max_count'}))
    mc_map = dict(zip(max_count['fusion_name'].str.upper(), max_count['mitelman_max_count']))
    mc_map.update(dict(zip(
        (max_count['fusion_name'].str.split('::').str[1]+'::'+
         max_count['fusion_name'].str.split('::').str[0]).str.upper(),
        max_count['mitelman_max_count']
    )))

    mit_by_cancer = {}
    for ct in mit['tcga_cancer_type'].unique():
        sub = mit[mit['tcga_cancer_type']==ct]
        sp  = sub['fusion_name'].str.split('::')
        m   = dict(zip(sub['fusion_name'].str.upper(), sub['Count']))
        m.update(dict(zip((sp.str[1]+'::'+sp.str[0]).str.upper(), sub['Count'])))
        mit_by_cancer[ct] = m

    fn_up = df['fusion_name'].str.upper()
    df['in_mitelman_exact']       = fn_up.isin(mit_pairs).astype(int)
    df['in_mitelman_same_cancer'] = df.apply(
        lambda r: int(r['cancer_type'] in mit_by_cancer and
                      r['fusion_name'].upper() in mit_by_cancer[r['cancer_type']]), axis=1)
    df['mitelman_literature_count']  = fn_up.map(mc_map).fillna(0).astype(int)
    df['mitelman_same_cancer_count'] = df.apply(
        lambda r: mit_by_cancer.get(r['cancer_type'],{}).get(r['fusion_name'].upper(),0), axis=1)
    df['mitelman_match_score'] = df.apply(
        lambda r: (2 if r['cancer_type'] in mit_by_cancer and
                        r['fusion_name'].upper() in mit_by_cancer[r['cancer_type']]
                   else 1 if r['fusion_name'].upper() in mit_pairs else 0), axis=1)

    print(f'    in_mitelman_exact: {df["in_mitelman_exact"].sum()}')
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [4]  CANCERMINE + HUGO GENE FUNCTION FLAGS                  ║
# ║       → gene1/2_{oncogene, tsg, driver, is_kinase, is_tf,    ║
# ║          is_rna_binding, is_chromatin, is_rtk,               ║
# ║          is_phosphatase, is_cell_cycle, is_ubiquitin,        ║
# ║          is_gtpase, is_dna_repair}                           ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_gene_function(df, refs):
    """
    Adds CancerMine roles (oncogene/tsg/driver) and HUGO function flags
    (kinase, tf, rna_binding, chromatin, rtk, phosphatase, cell_cycle,
     ubiquitin, gtpase, dna_repair).
    Used by: relevance score (tiers 1-3) + model features
    """
    print('CancerMine + HUGO gene function...')
    df = df.copy()

    # ── CancerMine ────────────────────────────────────────────
    cm = pd.read_csv(refs['cancermine_file'], sep='\t')
    cm['role'] = cm['role'].str.upper().str.strip()
    cm['gene_normalized'] = cm['gene_normalized'].str.upper().str.strip()

    roles = cm.groupby(['gene_normalized','role']).size().unstack(fill_value=0)
    roles = (roles > 0).astype(int)
    roles['DRIVER'] = 1

    onco_map   = roles['ONCOGENE'].to_dict()       if 'ONCOGENE'          in roles.columns else {}
    tsg_map    = roles['TUMOR_SUPPRESSOR'].to_dict() if 'TUMOR_SUPPRESSOR' in roles.columns else {}
    driver_map = roles['DRIVER'].to_dict()

    for g in [1, 2]:
        gs = df[f'gene{g}'].str.upper().str.strip()
        df[f'gene{g}_oncogene'] = gs.map(onco_map).fillna(0).astype(int)
        df[f'gene{g}_tsg']      = gs.map(tsg_map).fillna(0).astype(int)
        df[f'gene{g}_driver']   = gs.map(driver_map).fillna(0).astype(int)
        # keep legacy _is_oncogene alias used by relevance score
        df[f'gene{g}_is_oncogene'] = df[f'gene{g}_oncogene']
        df[f'gene{g}_is_tsg']      = df[f'gene{g}_tsg']

    # ── HUGO gene sets ────────────────────────────────────────
    with open(refs['hugo_gene_sets'], 'r') as f:
        gts = json.load(f)

    target_sets = {
        'kinase':      set(gts['kinase']),
        'tf':          set(gts['tf']),
        'rna_binding': set(gts['rna_binding']),
        'chromatin':   set(gts['chromatin']),
        'rtk':         set(gts['rtk']),
        'phosphatase': set(gts['phosphatase']),
        'cell_cycle':  set(gts['cell_cycle']),
        'ubiquitin':   set(gts['ubiquitin']),
        'gtpase':      set(gts['gtpase']),
        'dna_repair':  set(gts['dna_repair']),
    }

    new_cols = {}
    for g in [1, 2]:
        gs = df[f'gene{g}'].str.upper().str.strip()
        for flag, gene_set in target_sets.items():
            new_cols[f'gene{g}_is_{flag}'] = gs.isin(gene_set).astype(int)
    df = pd.concat([df, pd.DataFrame(new_cols, index=df.index)], axis=1)

    df['n_oncogenic_partners'] = df['gene1_oncogene'] + df['gene2_oncogene']
    df['has_tf']               = ((df['gene1_is_tf']==1) | (df['gene2_is_tf']==1)).astype(int)
    
    print(f'    oncogene g1/g2: {df["gene1_oncogene"].sum()} / {df["gene2_oncogene"].sum()} | '
          f'kinase: {(df["gene1_is_kinase"]+df["gene2_is_kinase"]>0).sum()}')
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [5]  CLINICAL PANELS                                        ║
# ║       → gene1/2_in_any_panel                                 ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_panels(df, refs):
    """
    Adds: gene1_in_any_panel, gene2_in_any_panel
    Model feature.
    """
    print('Clinical panels...')
    df = df.copy()
    panel_df = pd.read_csv(refs['fusion_panels'], sep='\t')
    any_panel = set(panel_df['gene'].str.upper())
    for g in [1, 2]:
        df[f'gene{g}_in_any_panel'] = df[f'gene{g}'].str.upper().str.strip().isin(any_panel).astype(int)
    print(f'    gene1_in_any_panel: {df["gene1_in_any_panel"].sum()} | '
          f'gene2_in_any_panel: {df["gene2_in_any_panel"].sum()}')
    
    df['panel_gene_count'] = df['gene1_in_any_panel'] + df['gene2_in_any_panel']
    
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [6]  TAD SAME DOMAIN                                        ║
# ║       → tad_same_domain                                      ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_tad(df, refs):
    """
    Adds: tad_same_domain
    Model feature.
    """
    print('TAD boundary...')
    df = df.copy()
    with open(refs['tad_pickle'], 'rb') as f:
        tad_dict = pickle.load(f)
    g2t = {g.strip().upper(): t for t, genes in tad_dict.items() for g in genes}
    t1 = df['gene1'].str.strip().map(g2t)
    t2 = df['gene2'].str.strip().map(g2t)
    df['tad_same_domain'] = ((t1 == t2) & t1.notna()).astype(int)
    print(f'    tad_same_domain=1: {df["tad_same_domain"].sum():,}')
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [7]  PPI NETWORK  (OmniPath)                                ║
# ║       → gene1/2_centrality, centrality_sum/max_log           ║
# ║       → ppi_direct_interaction, ppi_interaction_score        ║
# ║       → ppi_shared_neighbors, ppi_jaccard_similarity         ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_ppi(df, refs):
    """
    Adds: gene1/2_centrality, centrality_sum_log, centrality_max_log,
          ppi_direct_interaction, ppi_interaction_score,
          ppi_shared_neighbors, ppi_jaccard_similarity

    Reads pre-built cache files from explicit refs:
      refs['ppi_centrality_file'] — gene -> degree_centrality
      refs['ppi_edges_file']      — gene1, gene2, interaction_score
      refs['ppi_adjacency_file']  — {gene: set(neighbors)}

    Run preprocess_omnipath.py once to generate these files.
    """
    print('PPI network (from cache)...')
    df = df.copy()

    # ── Load cache files ──────────────────────────────────────
    cent_df = pd.read_csv(refs['ppi_centrality_file'], sep='\t')
    edge_df = pd.read_csv(refs['ppi_edges_file'], sep='\t')
    with open(refs['ppi_adjacency_file'], 'rb') as f:
        adj = pickle.load(f)

    print(f'    Loaded: {len(cent_df):,} genes | {len(edge_df):,} edges')

    # ── Build fast lookup dicts ───────────────────────────────
    cent_map = dict(zip(cent_df['gene'].str.upper(), cent_df['centrality']))

    # Edge set for direct interaction check + score lookup
    # Stored as gene1 < gene2 (sorted), so normalize both ways
    iset = set()
    sc   = {}
    for _, r in edge_df.iterrows():
        g1, g2 = r['gene1'], r['gene2']
        iset.add((g1, g2)); iset.add((g2, g1))
        sc[(g1, g2)] = sc[(g2, g1)] = r['interaction_score']

    # ── Annotate ──────────────────────────────────────────────
    g1u = df['gene1'].str.strip().str.upper()
    g2u = df['gene2'].str.strip().str.upper()

    s1 = g1u.map(cent_map).fillna(0.0)
    s2 = g2u.map(cent_map).fillna(0.0)

    df['gene1_centrality']   = s1
    df['gene2_centrality']   = s2
    df['centrality_sum_log'] = np.log1p(s1 + s2)
    df['centrality_max_log'] = np.log1p(np.maximum(s1, s2))

    pairs = list(zip(g1u, g2u))
    df['ppi_direct_interaction'] = [int(p in iset) for p in pairs]
    df['ppi_interaction_score']  = [sc.get(p, 0.0) for p in pairs]

    def _shared(a, b):
        n1, n2 = adj.get(a, set()), adj.get(b, set())
        sh = len(n1 & n2)
        u  = len(n1 | n2)
        return sh, sh / u if u else 0.0

    results = [_shared(a, b) for a, b in pairs]
    df['ppi_shared_neighbors']   = [r[0] for r in results]
    df['ppi_jaccard_similarity'] = [r[1] for r in results]

    print(f'    direct interactions: {df["ppi_direct_interaction"].sum():,}')
    return df


# ╔══════════════════════════════════════════════════════════════════╗
# ║  [8] DRUGGABILITY                                                ║
# ║       → fusion_is_druggable, fusion_is_druggable_cancer_specific ║
# ║       → fusion_has_druggable_gene, gene1/2_is_druggable          ║
# ║       → gene1/2_druggable_cancer_specific                        ║
# ║       → fusion_gene_druggable_cancer_specific                    ║
# ╚══════════════════════════════════════════════════════════════════╝

def annotate_druggability(df, refs):
    """
    Adds: fusion_is_druggable, fusion_is_druggable_cancer_specific,
          gene1/2_is_druggable, fusion_has_druggable_gene,
          gene1/2_druggable_cancer_specific,
          fusion_gene_druggable_cancer_specific
    Used by: relevance score (all tiers)
    """
    print('Druggability...')
    df = df.copy()
    drug = pd.read_csv(refs['druggable_file'], sep='\t')

    fusion_cancer_map = {}
    for _, r in drug.iterrows():
        g1 = str(r['gene1']).strip().upper(); g2 = str(r['gene2']).strip().upper()
        ct = str(r['all_actionable_cancers']).strip().upper()
        for k in [(g1,g2),(g2,g1)]:
            fusion_cancer_map.setdefault(k, set()).add(ct)

    gene_cancer_map = {}
    for _, r in drug.iterrows():
        g  = str(r['gene_targets']).strip().upper()
        ct = str(r['all_actionable_cancers']).strip().upper()
        gene_cancer_map.setdefault(g, set()).add(ct)

    def _check(row):
        g1 = str(row['gene1']).strip().upper(); g2 = str(row['gene2']).strip().upper()
        ct = str(row['cancer_type']).strip().upper()
        fk = (g1,g2) if (g1,g2) in fusion_cancer_map else (g2,g1) if (g2,g1) in fusion_cancer_map else None
        fid  = 1 if fk else 0
        fics = 1 if fk and ct in (fusion_cancer_map[fk] - {'ALLTYPES'}) else 0
        g1d  = 1 if g1 in gene_cancer_map else 0
        g2d  = 1 if g2 in gene_cancer_map else 0
        fhd  = 1 if g1d or g2d else 0
        g1cs = 1 if g1 in gene_cancer_map and ct in (gene_cancer_map[g1] - {'ALLTYPES'}) else 0
        g2cs = 1 if g2 in gene_cancer_map and ct in (gene_cancer_map[g2] - {'ALLTYPES'}) else 0
        fgcs = 1 if g1cs or g2cs else 0
        return pd.Series([fid,fics,g1d,g2d,fhd,g1cs,g2cs,fgcs])

    cols = ['fusion_is_druggable','fusion_is_druggable_cancer_specific',
            'gene1_is_druggable','gene2_is_druggable','fusion_has_druggable_gene',
            'gene1_druggable_cancer_specific','gene2_druggable_cancer_specific',
            'fusion_gene_druggable_cancer_specific']
    df[cols] = df.apply(_check, axis=1)
    print(f'    fusion_is_druggable: {df["fusion_is_druggable"].sum()} | '
          f'cancer_specific: {df["fusion_is_druggable_cancer_specific"].sum()}')
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [9] BREAKPOINT REGION ANNOTATION  (bedtools)                ║
# ║       → gene1/2_region_type, gene1/2_bkp_category            ║
# ║       → gene1/2_transcript_id  (needed by Part 12)           ║
# ╚══════════════════════════════════════════════════════════════╝

_BKP_COLS = ['bp_chr','bp_start','bp_end','bp_name','bp_score','bp_strand',
             'fusion_name','gene1_symbol','gene2_symbol',
             'gene1_chr','gene1_bkp','gene2_chr','gene2_bkp',
             'gene1_strand','gene2_strand',
             'ann_chr','ann_start','ann_end','ann_region_type',
             'ann_gene_name','ann_strand','ann_gene_id',
             'ann_transcript_id','ann_exon_number','ann_mane_select']

_LOCI_COLS = ['bp_chr','bp_start','bp_end','which_gene','bp_score','bp_strand',
              'fusion_name','gene1_symbol','gene2_symbol',
              'gene1_chr','gene1_bkp','gene2_chr','gene2_bkp',
              'gene1_strand','gene2_strand',
              'ann_chr','ann_start','ann_end','ann_gene_name','ann_score','ann_strand']


_BKP_OUTPUT_COLS = ['gene1_region_type', 'gene1_bkp_category', 'gene1_transcript_id',
                    'gene2_region_type', 'gene2_bkp_category', 'gene2_transcript_id']


def annotate_breakpoints(df, refs, output_dir, project_name):
    """
    Adds: gene1/2_region_type, gene1/2_bkp_category, gene1/2_transcript_id
    Also corrects gene names via GENCODE loci overlap.
    Model features. Part 12 (cDNA) depends on this running first.

    If breakpoint columns (gene1/2_chr, gene1/2_bkp, gene1/2_strand) are absent
    from df, all output columns are set to pd.NA and bedtools is skipped entirely.
    """
    print('Breakpoint annotation (bedtools)...')

    if not _has_breakpoint_cols(df):
        print('    ⚠️  Breakpoint columns not found — skipping bedtools, filling outputs with NA')
        df = df.copy()
        for col in _BKP_OUTPUT_COLS:
            df[col] = pd.NA
        return df

    df = df.copy()
    os.makedirs(output_dir, exist_ok=True)

    # Write BED (only needs the core positional columns)
    bed_path = os.path.join(output_dir, f'{project_name}_breakpoints.bed')
    rows = []
    for side, cc, bc, sc in [('gene1','gene1_chr','gene1_bkp','gene1_strand'),
                               ('gene2','gene2_chr','gene2_bkp','gene2_strand')]:
        rows.append(pd.DataFrame({
            'chr':df[cc], 'start':df[bc]-1, 'end':df[bc],
            'name':side, 'score':'.', 'strand':df[sc],
            'fusion_name':df['fusion_name'],
            'gene1_symbol':df['gene1'], 'gene2_symbol':df['gene2'],
            'gene1_chr':df['gene1_chr'], 'gene1_bkp':df['gene1_bkp'],
            'gene2_chr':df['gene2_chr'], 'gene2_bkp':df['gene2_bkp'],
            'gene1_strand':df['gene1_strand'], 'gene2_strand':df['gene2_strand'],
        }))
    bed_df = pd.concat(rows).copy()
    bed_df['chr'] = bed_df['chr'].astype(str).str.strip()
    bed_df['start'] = pd.to_numeric(bed_df['start'], errors='coerce')
    bed_df['end'] = pd.to_numeric(bed_df['end'], errors='coerce')

    # Guard against malformed coordinates that make bedtools fail.
    n_before = len(bed_df)
    invalid = (
        bed_df['chr'].eq('') |
        bed_df['start'].isna() |
        bed_df['end'].isna() |
        (bed_df['end'] <= 0) |
        (~bed_df['strand'].isin(['+', '-']))
    )
    if invalid.any():
        print(f"    Dropping {int(invalid.sum()):,} malformed breakpoint rows before bedtools")
    bed_df = bed_df[~invalid].copy()
    bed_df['end'] = bed_df['end'].astype(int)
    bed_df['start'] = bed_df['start'].astype(int).clip(lower=0)
    bed_df = bed_df[bed_df['end'] > bed_df['start']].copy()
    if len(bed_df) < n_before:
        print(f"    BED rows kept: {len(bed_df):,}/{n_before:,}")

    bed_df = bed_df.sort_values(['chr', 'start'])
    bed_df.to_csv(bed_path, sep='	', header=False, index=False)

    a_bt = pybedtools.BedTool.from_dataframe(bed_df)

    def _run(b_ref, out):
        try:
            b_bt = pybedtools.BedTool(b_ref)
            result = a_bt.intersect(b_bt, s=True, wa=True, wb=True, loj=True)
            result.saveas(out)
        except Exception as exc:
            raise RuntimeError(
                f"pybedtools intersect failed for reference: {b_ref}\n"
                f"error: {exc}"
            ) from exc

    # ── Loci overlap -> correct gene names ────────────────────
    # Work on a slim positional df so merges don't touch the full df
    core_cols = ['fusion_name','gene1','gene2',
                 'gene1_chr','gene1_bkp','gene1_strand',
                 'gene2_chr','gene2_bkp','gene2_strand']
    pos = df[core_cols].copy()

    loci_out = os.path.join(output_dir, f'{project_name}_loci.bed')
    _run(refs['gencode_loci_bed'], loci_out)
    dl = pd.read_csv(loci_out, sep='\t', header=None, names=_LOCI_COLS)
    dl['bp_pos'] = dl['bp_end']
    dl = dl.drop_duplicates(subset=['bp_chr','bp_pos','bp_strand','which_gene'])

    for side in ['gene1','gene2']:
        sub = dl[dl['which_gene']==side][
            ['bp_chr','bp_pos','bp_strand','ann_gene_name']
        ].rename(columns={
            'bp_chr':f'{side}_chr', 'bp_pos':f'{side}_bkp',
            'bp_strand':f'{side}_strand', 'ann_gene_name':f'{side}_gencode'
        })
        pos = pos.merge(sub, on=[f'{side}_chr',f'{side}_bkp',f'{side}_strand'], how='left')
        pos[f'{side}_gencode'] = pos[f'{side}_gencode'].apply(
            lambda x: None if pd.isna(x) or str(x).startswith('ENSG') or str(x)=='.' else x)
        # Update gene name in pos only; write corrected names back to df
        pos[side] = pos[f'{side}_gencode'].fillna(pos[side])
        pos = pos.drop(columns=[f'{side}_gencode'])

    # Write corrected gene names back to the full df
    df['gene1'] = pos['gene1'].values
    df['gene2'] = pos['gene2'].values

    # ── Detailed BED -> region type + transcript ──────────────
    bkp_out = os.path.join(output_dir, f'{project_name}_bkp_annot.bed')
    _run(refs['gencode_bkp_bed'], bkp_out)
    db = pd.read_csv(bkp_out, sep='\t', header=None, names=_BKP_COLS)
    for c in ['ann_chr','ann_start','ann_end','ann_region_type','ann_gene_name',
              'ann_strand','ann_gene_id','ann_transcript_id','ann_exon_number','ann_mane_select']:
        db[c] = db[c].replace(['.', '-1'], pd.NA)
    for c in ['bp_start','bp_end','gene1_bkp','gene2_bkp']:
        db[c] = pd.to_numeric(db[c], errors='coerce').astype('Int64')
    db['ann_exon_number'] = pd.to_numeric(db['ann_exon_number'], errors='coerce').astype('Int64')
    db['ann_mane_select'] = (db['ann_mane_select'].astype(str).str.lower()
                              .map({'true':True,'false':False}))

    # Remove strand/gene mismatches
    ann = db[db['ann_chr'].notna()].copy()
    g1b = ((ann['bp_name']=='gene1') &
           ((ann['ann_gene_name']!=ann['gene1_symbol']) |
            (ann['ann_strand']   !=ann['gene1_strand'])))
    g2b = ((ann['bp_name']=='gene2') &
           ((ann['ann_gene_name']!=ann['gene2_symbol']) |
            (ann['ann_strand']   !=ann['gene2_strand'])))
    db  = db[~db.index.isin(ann[g1b|g2b].index)]

    # Prioritize: MANE > CDS > other
    def _pick(group):
        if group['ann_chr'].isna().all(): return group.iloc[0]
        a    = group[group['ann_chr'].notna()]
        mane = a[a['ann_mane_select']==True]
        if not mane.empty: return mane.iloc[0]
        cds  = a[a['ann_region_type']=='CDS']
        if not cds.empty: return cds.iloc[0]
        return a.iloc[0]

    db_pri = (db.groupby(['fusion_name','gene1_bkp','gene2_bkp','bp_name'], group_keys=True)
               .apply(_pick, include_groups=False).reset_index())

    # Nearest exon boundary for introns
    gencode_bed = pd.read_csv(refs['gencode_bkp_bed'], sep='\t', header=None,
                               names=['chr','start','end','region_type','gene_name','strand',
                                      'gene_id','transcript_id','exon_number','mane_select'])
    exons = gencode_bed[gencode_bed['region_type'].isin(['CDS','UTR5','UTR3'])].copy()
    exons['s1'] = exons['start'].astype(int) + 1
    exons['e1'] = exons['end'].astype(int)

    def _nearest(bdf, side):
        bdf[f'{side}_nearest_dist'] = pd.NA
        has = bdf['ann_transcript_id'].notna()
        if not has.any(): return bdf
        sub = bdf[has].copy().reset_index().rename(columns={'index':'oi'})
        m   = sub.merge(exons[['transcript_id','s1','e1']],
                        left_on='ann_transcript_id', right_on='transcript_id',
                        how='left').dropna(subset=['s1'])
        if m.empty: return bdf
        bpos = m[f'{side}_bkp'].astype(int)
        m['dist'] = np.minimum((bpos - m['s1']).abs(), (bpos - m['e1']).abs())
        best = m.loc[m.groupby('oi')['dist'].idxmin(), ['oi','dist']]
        bdf.loc[best['oi'], f'{side}_nearest_dist'] = best['dist'].values
        return bdf

    def _reclass(rt, dist, thr=200):
        if pd.isna(rt): return 'Intergenic'
        rt = str(rt)
        if rt in ('CDS','UTR5','UTR3'): return rt
        if rt in ('gene','intergenic'):  return 'Intergenic'
        if rt == 'intron':
            if pd.isna(dist): return 'Intergenic'
            return 'Near_Exon_Intron' if float(dist) <= thr else 'Deep_Intron'
        return 'Intergenic'

    # Merge breakpoint annotations back onto the full df (don't reassign df to slim subset)
    for side in ['gene1','gene2']:
        sdf = db_pri[db_pri['bp_name']==side].copy()
        sdf = _nearest(sdf, side)
        sdf[f'{side}_region_type']   = sdf['ann_region_type']
        sdf[f'{side}_bkp_category']  = sdf.apply(
            lambda r: _reclass(r['ann_region_type'], r.get(f'{side}_nearest_dist', pd.NA)), axis=1)
        sdf[f'{side}_transcript_id'] = sdf['ann_transcript_id']
        mc = ['fusion_name','gene1_bkp','gene2_bkp',
              f'{side}_region_type', f'{side}_bkp_category', f'{side}_transcript_id']
        df = df.merge(sdf[mc], on=['fusion_name','gene1_bkp','gene2_bkp'], how='left')

    print(f'    gene1_bkp_category: {df["gene1_bkp_category"].value_counts().to_dict()}')
    return df



# ╔══════════════════════════════════════════════════════════════╗
# ║  [10] cDNA / PYENSEMBL FEATURES                              ║
# ║       → gene1/2_frame_mod3, fusion_in_frame, has_early_stop  ║
# ║       → gene1/2_fusion_pct                                   ║
# ║       → log_gene1/2_len_nt, log_gene1/2_utr5/3_len           ║
# ╚══════════════════════════════════════════════════════════════╝

_CDNA_OUTPUT_COLS = [
    'gene1_frame_mod3', 'gene2_frame_mod3',
    'fusion_in_frame', 'has_early_stop',
    'gene1_fusion_pct', 'gene2_fusion_pct',
    'log_gene1_len_nt', 'log_gene2_len_nt',
    'log_gene1_utr5_len', 'log_gene1_utr3_len',
    'log_gene2_utr5_len', 'log_gene2_utr3_len',
]

def annotate_cdna_features(df, refs):
    """
    Requires: gene1/2_transcript_id  (from Part 11)
    Adds: gene1/2_frame_mod3, fusion_in_frame, has_early_stop,
          gene1/2_fusion_pct, log_gene1/2_len_nt,
          log_gene1/2_utr5_len, log_gene1/2_utr3_len
    Model features.

    If gene1/2_transcript_id columns are absent or entirely NA (i.e. [11]
    was skipped due to missing breakpoint columns), all output columns are
    set to pd.NA and the function returns early.
    """
    print('cDNA features (PyEnsembl)...')

    # Guard: transcript_id columns must exist and have at least one non-NA value
    tx_cols_present = (
        'gene1_transcript_id' in df.columns and
        'gene2_transcript_id' in df.columns
    )
    tx_cols_populated = tx_cols_present and (
        df['gene1_transcript_id'].notna().any() or
        df['gene2_transcript_id'].notna().any()
    )

    if not tx_cols_populated:
        print('    ⚠️  gene1/2_transcript_id not available — skipping cDNA annotation, filling outputs with NA')
        df = df.copy()
        for col in _CDNA_OUTPUT_COLS:
            df[col] = pd.NA
        return df

    import pyensembl
    from Bio.Seq import Seq
    df = df.copy()
    gd = pyensembl.Genome(
        reference_name='GRCh38', annotation_name='gencode', annotation_version='v49',
        gtf_path_or_url=refs['gencode_gtf'],
        transcript_fasta_paths_or_urls=refs['gencode_tx_fasta'],
        protein_fasta_paths_or_urls=refs['gencode_prot_fasta'],
        cache_directory_path=refs['gencode_cache_dir'],
        copy_local_files_to_cache=False,
    )
    gd.index()
    tx_ids = list(set(df['gene1_transcript_id'].dropna().tolist() +
                      df['gene2_transcript_id'].dropna().tolist()))
    tx_cache = {}
    for tx in tx_ids:
        try: tx_cache[tx] = gd.transcript_by_id(tx)
        except: tx_cache[tx] = None
    print(f'    Cached {sum(v is not None for v in tx_cache.values()):,}/{len(tx_cache):,} transcripts')
    txinfo = {}
    for tx_id, tx in tx_cache.items():
        if tx is None: txinfo[tx_id] = None; continue
        smap = {}; offset = 1; minus = tx.strand == '-'
        exons = sorted(tx.exons, key=lambda e: e.start, reverse=minus)
        for ex in exons:
            rng = range(ex.end, ex.start-1, -1) if minus else range(ex.start, ex.end+1)
            for pos in rng: smap[pos] = offset; offset += 1
        try: utr5 = len(tx.five_prime_utr_sequence)  if tx.five_prime_utr_sequence  else 0
        except: utr5 = 0
        try: utr3 = len(tx.three_prime_utr_sequence) if tx.three_prime_utr_sequence else 0
        except: utr3 = 0
        try: coding = str(tx.coding_sequence) if tx.coding_sequence else ''
        except: coding = ''
        txinfo[tx_id] = dict(smap=smap, utr5=utr5, utr3=utr3, coding=coding)
    def _extract(row):
        def _part(tx_id, bkp, is_g1, bkp_category=None):
            info = txinfo.get(tx_id)
            if not info: return '', None
            bkp_int = int(bkp)
            off = info['smap'].get(bkp_int)
            # If breakpoint falls in intron, snap to nearest exon boundary
            if off is None and bkp_category in ('Near_Exon_Intron',):
                nearest = min(info['smap'].keys(), key=lambda p: abs(p - bkp_int))
                off = info['smap'].get(nearest)
            if off is None: return '', None
            coff = off - info['utr5']
            seq = info['coding'][:coff] if is_g1 else info['coding'][coff-1:]
            return seq.upper(), len(seq) % 3
        g1s, g1m = _part(row.get('gene1_transcript_id'), row['gene1_bkp'], True,  row.get('gene1_bkp_category'))
        g2s, g2m = _part(row.get('gene2_transcript_id'), row['gene2_bkp'], False, row.get('gene2_bkp_category'))
        fused = g1s + g2s
        in_frame = (1 if len(fused) % 3 == 0 else 0) if (g1s and g2s) else 1
        hes = pd.NA
        if fused:
            prot = str(Seq(fused[:len(fused)-len(fused)%3]).translate())
            fs = prot.find('*')
            hes = int(fs != -1 and fs < len(prot)-1)
        return pd.Series([len(g1s), len(g2s), g1m, g2m, in_frame, hes])
    df[['gene1_len_nt','gene2_len_nt','gene1_frame_mod3','gene2_frame_mod3',
        'fusion_in_frame','has_early_stop']] = df.apply(_extract, axis=1)
    def _pct(row):
        res = [0.0, 0.0]
        for i, side in enumerate(['gene1','gene2']):
            info = txinfo.get(row.get(f'{side}_transcript_id'))
            if info and info['coding']:
                res[i] = row[f'{side}_len_nt'] / len(info['coding']) * 100
        return pd.Series(res)
    df[['gene1_fusion_pct','gene2_fusion_pct']] = df.apply(_pct, axis=1)
    def _utr(row):
        vals = []
        for side in ['gene1','gene2']:
            info = txinfo.get(row.get(f'{side}_transcript_id'))
            vals += [info['utr5'], info['utr3']] if info else [0, 0]
        return pd.Series(vals)
    df[['gene1_utr5_len','gene1_utr3_len','gene2_utr5_len','gene2_utr3_len']] = df.apply(_utr, axis=1)
    df['log_gene1_len_nt']   = np.log1p(df['gene1_len_nt'])
    df['log_gene2_len_nt']   = np.log1p(df['gene2_len_nt'])
    df['log_gene1_utr5_len'] = np.log1p(df['gene1_utr5_len'])
    df['log_gene1_utr3_len'] = np.log1p(df['gene1_utr3_len'])
    df['log_gene2_utr5_len'] = np.log1p(df['gene2_utr5_len'])
    df['log_gene2_utr3_len'] = np.log1p(df['gene2_utr3_len'])
    in_frame_pct = (df['fusion_in_frame']=='Yes').mean()*100
    print(f'    in_frame: {in_frame_pct:.1f}% | early_stop: {df["has_early_stop"].sum()}')

    df = df.drop(columns=['gene1_len_nt', 'gene2_len_nt',
                          'gene1_utr5_len', 'gene1_utr3_len', 'gene2_utr5_len', 'gene2_utr3_len'])
    return df


# ╔══════════════════════════════════════════════════════════════╗
# ║  [11] CANCER DB                                              ║
# ║       → gene1/2_in_cancer_db, gene1/2_cancer_db_mismatch     ║
# ║       → gene1/2_unknown_in_cancer                            ║
# ╚══════════════════════════════════════════════════════════════╝

# Valid TCGA + pediatric cancer type codes used to filter all three DBs
VALID_CANCER_TYPES = [
    # Hematologic
    'ALL', 'LAML', 'DLBC', 'UVM',
    # Solid tumors
    'BRCA', 'STAD', 'OV', 'LUSC', 'HNSC', 'ESCA', 'LGG',
    'LUAD', 'PRAD', 'KIRC', 'SKCM', 'THCA', 'SARC', 'BLCA',
    'GBM', 'LIHC', 'CESC', 'UCEC', 'KIRP', 'COAD', 'PCPG',
    'PAAD', 'KICH', 'READ', 'THYM', 'UCS', 'TGCT', 'ACC',
    'MESO', 'CHOL',
    # Pediatric
    'NBL', 'WT', 'OS', 'RT', 'CCSK',
]


def _build_cancer_db_lookup(refs, valid_cancer_types=VALID_CANCER_TYPES):
    """
    Build a unified gene → set(cancer_types) lookup using a tiered strategy:

    Tier A — COSMIC genes: use COSMIC cancer types only (high specificity).
    Tier B — Non-COSMIC genes from druggable fusions only (CIViC excluded
              due to pan-cancer noise from therapeutic evidence entries).
              Also filters out 'alltypes' entries.

    Returns: dict {GENE_UPPER: set of cancer type strings}
    """
    cosmic_df = pd.read_csv(refs['cosmic_file'], sep='	')
    drug_df   = pd.read_csv(refs['druggable_file'], sep='	')

    # ── Tier A: COSMIC ────────────────────────────────────────
    cosmic_pairs = pd.concat([
        cosmic_df[['fusion_name','tcga_cancer_type']].assign(
            gene=cosmic_df['fusion_name'].str.split('::').str[0].str.upper()),
        cosmic_df[['fusion_name','tcga_cancer_type']].assign(
            gene=cosmic_df['fusion_name'].str.split('::').str[1].str.upper()),
    ])[['gene','tcga_cancer_type']].rename(columns={'tcga_cancer_type':'cancer_type'})
    cosmic_pairs  = cosmic_pairs[cosmic_pairs['cancer_type'].isin(valid_cancer_types)]
    cosmic_genes  = set(cosmic_pairs['gene'].unique())
    cosmic_lookup = cosmic_pairs.groupby('gene')['cancer_type'].apply(set).to_dict()

    # ── Tier B: Druggable fusions only (CIViC removed — too noisy) ───
    drug_pairs = pd.concat([
        drug_df[['gene1','all_actionable_cancers']].rename(
            columns={'gene1':'gene','all_actionable_cancers':'cancer_type'}
        ).assign(gene=drug_df['gene1'].str.upper()),
        drug_df[['gene2','all_actionable_cancers']].rename(
            columns={'gene2':'gene','all_actionable_cancers':'cancer_type'}
        ).assign(gene=drug_df['gene2'].str.upper()),
    ])
    # Filter out 'alltypes' and invalid cancer types
    drug_pairs = drug_pairs[
        drug_pairs['cancer_type'].isin(valid_cancer_types) &
        (drug_pairs['cancer_type'].str.upper() != 'ALLTYPES')
    ]
    drug_lookup = drug_pairs.groupby('gene')['cancer_type'].apply(set).to_dict()

    # ── Combine ───────────────────────────────────────────────
    non_cosmic_genes = set(drug_lookup) - cosmic_genes
    unified = {}
    for gene, cancers in cosmic_lookup.items():   # Tier A
        unified[gene] = cancers
    for gene in non_cosmic_genes:                 # Tier B
        cancers = {ct for ct in drug_lookup.get(gene, set()) if ct.upper() != 'ALLTYPES'}
        if cancers:
            unified[gene] = cancers

    print(f'    COSMIC genes (Tier A): {len(cosmic_genes):,}')
    print(f'    Non-COSMIC genes (Tier B): {len(non_cosmic_genes):,}')
    print(f'    Total genes in lookup: {len(unified):,}')
    return unified


def annotate_cancer_db(df, refs, valid_cancer_types=VALID_CANCER_TYPES):
    """
    Builds a tiered gene→cancer-type lookup from COSMIC (Tier A),
    druggable fusions, and CIViC (Tier B), then annotates each fusion.

    Adds:
      gene1/2_in_cancer_db       — gene known in THIS cancer type
      gene1/2_cancer_db_mismatch — gene known, but NOT in this cancer type
      gene1/2_unknown_in_cancer  — gene not found in any cancer DB
    Model features.
    """
    print('Cancer DB annotation...')
    df = df.copy()

    gene_cancer_lookup = _build_cancer_db_lookup(refs, valid_cancer_types)
    all_known = set(gene_cancer_lookup)

    for g in [1, 2]:
        g_col = f'gene{g}'
        df[f'gene{g}_in_cancer_db'] = df.apply(
            lambda r: int(r['cancer_type'] in gene_cancer_lookup.get(r[g_col].upper(), set())),
            axis=1)
        df[f'gene{g}_cancer_db_mismatch'] = df.apply(
            lambda r: int(r[g_col].upper() in gene_cancer_lookup and
                          r['cancer_type'] not in gene_cancer_lookup.get(r[g_col].upper(), set())),
            axis=1)
        df[f'gene{g}_unknown_in_cancer'] = (~df[g_col].str.upper().isin(all_known)).astype(int)

    df['n_cancer_db_partners'] = df['gene1_in_cancer_db'] + df['gene2_in_cancer_db']
    
    for col in ['gene1_in_cancer_db','gene2_in_cancer_db',
                'gene1_cancer_db_mismatch','gene2_cancer_db_mismatch',
                'gene1_unknown_in_cancer','gene2_unknown_in_cancer']:
        print(f'    {col}: {df[col].sum():,} ({df[col].mean()*100:.1f}%)')
    return df



# ╔══════════════════════════════════════════════════════════════╗
# ║  [12] EXPRESSION  ← STUB                                     ║
# ║       → g1/2_log_exp, gene1/2_exp_low                        ║
# ║       → gene1/2_low_everywhere                               ║
# ║       → fusion_tissue_fit, exp_cancer_mismatch               ║
# ╚══════════════════════════════════════════════════════════════╝

def annotate_expression(df, refs, low_percentile=0.50, abs_floor=1.0):
    """
    Annotate fusions with HPA cancer-type specific expression features.
    Reads pre-built cache files from refs['hpa_cache_dir']:
      hpa_exp_dict.pkl       - {(GENE, CANCER_TYPE): median_ptpm}
      hpa_thresholds.tsv     - cancer_type, low_threshold
      hpa_low_everywhere.txt - genes low in all cancer types

    Run preprocess_hpa_expression.py once to generate these files.
    Falls back to computing from refs['hpa_file'] if cache not found.
    """
    print(f'Expression annotation (HPA, p{low_percentile*100:.0f} + floor={abs_floor} pTPM)...')
    df = df.copy()

    cache_dir = refs.get('hpa_cache_dir')
    exp_dict_path = refs.get('hpa_exp_dict_file')
    thresh_path = refs.get('hpa_thresholds_file')
    low_evry_path = refs.get('hpa_low_everywhere_file')

    # Backward-compatible fallback: derive cache paths from hpa_cache_dir.
    if cache_dir:
        exp_dict_path = exp_dict_path or os.path.join(cache_dir, 'hpa_exp_dict.pkl')
        thresh_path = thresh_path or os.path.join(cache_dir, 'hpa_thresholds.tsv')
        low_evry_path = low_evry_path or os.path.join(cache_dir, 'hpa_low_everywhere.txt')

    if all(p and os.path.exists(p) for p in [exp_dict_path, thresh_path, low_evry_path]):
        # Load from cache
        print('    Loading from cache...')
        with open(exp_dict_path, 'rb') as f:
            exp_dict = pickle.load(f)
        cancer_thresholds = pd.read_csv(thresh_path, sep='\t')
        with open(low_evry_path, 'r') as f:
            low_everywhere = set(f.read().splitlines())
        print(f'    {len(exp_dict):,} gene-cancer pairs | '
              f'{len(cancer_thresholds):,} cancer types | '
              f'{len(low_everywhere):,} low-everywhere genes')
    else:
        # Fall back: compute on the fly (slow)
        print('    Cache not found -- computing from raw HPA file (run preprocess_hpa_expression.py to cache)...')
        hpa_file = refs.get('hpa_file')
        if not hpa_file:
            raise FileNotFoundError(
                'HPA cache files are missing and no raw HPA file was found. '
                'Place a raw HPA table (e.g., hpa_expression.tsv) in --refs-dir '
                'or precompute cache files in hpa_cache_dir.'
            )
        hpa_ref = pd.read_csv(hpa_file)
        hpa_ref['gene_symbol'] = hpa_ref['gene_symbol'].str.upper()
        hpa_ref = hpa_ref.dropna(subset=['gene_symbol', 'cancer_type', 'median_ptpm'])
        hpa_ref = hpa_ref.drop_duplicates(subset=['gene_symbol', 'cancer_type'])
        exp_dict = hpa_ref.set_index(['gene_symbol', 'cancer_type'])['median_ptpm'].to_dict()
        cancer_thresholds = (
            hpa_ref.groupby('cancer_type')['median_ptpm']
            .quantile(low_percentile).rename('low_threshold').reset_index()
        )
        cancer_thresholds['low_threshold'] = cancer_thresholds['low_threshold'].clip(lower=abs_floor)
        hpa_with_thresh = hpa_ref.merge(cancer_thresholds, on='cancer_type', how='left')
        gene_low_flags = (
            hpa_with_thresh.groupby('gene_symbol')
            .apply(lambda g: (g['median_ptpm'] < g['low_threshold']).all(),
                   include_groups=False)
        )
        low_everywhere = set(gene_low_flags[gene_low_flags].index)

    threshold_dict = dict(zip(cancer_thresholds['cancer_type'],
                               cancer_thresholds['low_threshold']))

    # Per-fusion lookup
    g1u = df['gene1'].str.upper()
    g2u = df['gene2'].str.upper()

    df['gene1_exp'] = [exp_dict.get((g, ct), np.nan) for g, ct in zip(g1u, df['cancer_type'])]
    df['gene2_exp'] = [exp_dict.get((g, ct), np.nan) for g, ct in zip(g2u, df['cancer_type'])]

    g1_miss = df['gene1_exp'].isna().sum()
    g2_miss = df['gene2_exp'].isna().sum()
    print(f'    gene1 missing: {g1_miss:,} ({g1_miss/len(df)*100:.1f}%) -- filled with 0')
    print(f'    gene2 missing: {g2_miss:,} ({g2_miss/len(df)*100:.1f}%) -- filled with 0')

    df['gene1_exp'] = df['gene1_exp'].fillna(0)
    df['gene2_exp'] = df['gene2_exp'].fillna(0)
    df['_low_threshold'] = df['cancer_type'].map(threshold_dict)

    # Derived features
    df['g1_log_exp']           = np.log1p(df['gene1_exp'])
    df['g2_log_exp']           = np.log1p(df['gene2_exp'])
    df['gene1_exp_low']        = (df['gene1_exp'] < df['_low_threshold']).astype(int)
    df['gene2_exp_low']        = (df['gene2_exp'] < df['_low_threshold']).astype(int)
    df['gene1_low_everywhere'] = g1u.isin(low_everywhere).astype(int)
    df['gene2_low_everywhere'] = g2u.isin(low_everywhere).astype(int)
    df['fusion_tissue_fit']    = (
        (df['gene1_exp'] >= df['_low_threshold']) &
        (df['gene2_exp'] >= df['_low_threshold'])
    ).astype(int)
    df['exp_cancer_mismatch']  = (
        (df['gene1_exp_low'] == 1) ^ (df['gene2_exp_low'] == 1)
    ).astype(int)

    df['n_expressed_partners'] = (
        (df['gene1_exp_low']==0).astype(int) + (df['gene2_exp_low']==0).astype(int)
    )
    
    df = df.drop(columns=['_low_threshold', 'gene1_exp', 'gene2_exp'])

    for col in ['fusion_tissue_fit', 'gene1_exp_low', 'gene2_exp_low',
                'gene1_low_everywhere', 'gene2_low_everywhere', 'exp_cancer_mismatch']:
        print(f'    {col}: {df[col].sum():,} ({df[col].mean()*100:.1f}%)')
    return df


def run_full_pipeline(input_file, output_dir, project_name, refs):
    """Run all annotation steps and save final TSV."""
    os.makedirs(output_dir, exist_ok=True)
    df = preprocess_fusion_df(input_file)
    df = annotate_civic(df, refs)
    df = annotate_cosmic(df, refs)
    df = annotate_mitelman(df, refs)
    df = annotate_gene_function(df, refs)
    df = annotate_panels(df, refs)
    df = annotate_tad(df, refs)
    df = annotate_ppi(df, refs)
    df = annotate_druggability(df, refs)
    df = annotate_breakpoints(df, refs, output_dir, project_name)
    df = annotate_cdna_features(df, refs)
    df = annotate_cancer_db(df, refs)
    df = annotate_expression(df, refs, low_percentile=0.50)

    out = os.path.join(output_dir, f'{project_name}_annotated.tsv')
    df.to_csv(out, sep='\t', index=False)
    print(f'\n✅  Saved {len(df):,} fusions → {out}')
    return df



# ╔══════════════════════════════════════════════════════════════╗
# ║  [13]  RANKING                                               ║
# ╚══════════════════════════════════════════════════════════════╝


def load_model(model_dir):
    import lightgbm as lgb
    import pickle
    model = lgb.Booster(model_file=f'{model_dir}/fusion_lambdarank.txt')
    with open(f'{model_dir}/feature_cols.pkl', 'rb') as f:
        feature_cols = pickle.load(f)
    with open(f'{model_dir}/categorical_cols.pkl', 'rb') as f:
        categorical_cols = pickle.load(f)
    print(f'✅ Model loaded: {len(feature_cols)} features')
    return model, feature_cols, categorical_cols


def rank_fusions(df, model, feature_cols, categorical_cols):
    working_df = df.copy()
    X = working_df[feature_cols].copy()

    # Align categorical dtypes exactly with model training metadata.
    model_categories = getattr(model, 'pandas_categorical', None) or []
    category_map = {
        col: cats
        for col, cats in zip(categorical_cols, model_categories)
        if col in X.columns
    }

    for col in X.columns:
        if col in category_map:
            X[col] = X[col].astype('string').astype('category').cat.set_categories(category_map[col])
        elif str(X[col].dtype) in {'object', 'category', 'string'}:
            # Non-categorical model features should be numeric. Coerce defensively.
            X[col] = pd.to_numeric(X[col], errors='coerce').fillna(0)

    working_df['predicted_score'] = model.predict(X)
    working_df['rank'] = (
        working_df.groupby('sample_id')['predicted_score']
        .rank(ascending=False, method='first')
        .astype(int)
    )
    return working_df.sort_values(['sample_id', 'rank']).reset_index(drop=True)


def make_ranked_input_tsv(df, input_file):
    """Return ranked output with original input columns, Rank, and predicted_score."""
    input_cols = pd.read_csv(input_file, sep='\t', nrows=0).columns.tolist()
    output_cols = [
        col for col in input_cols
        if col not in ['Rank', 'predicted_score']
    ] + ['Rank', 'predicted_score']

    ranked_df = df.copy()
    ranked_df['Rank'] = ranked_df['rank'].astype(int)

    missing_cols = [
        col for col in input_cols
        if col not in ranked_df.columns
        and col not in ['Rank', 'predicted_score']
    ]
    if missing_cols:
        raise ValueError(
            'Ranked output cannot be created because input columns are missing after processing: '
            + ', '.join(missing_cols)
        )

    return ranked_df.sort_values(
        ['sample_id', 'Rank']
    ).loc[:, output_cols].reset_index(drop=True)


# ╔══════════════════════════════════════════════════════════════╗
# ║  [14]  SUMMARY OUTPUT                                        ║
# ╚══════════════════════════════════════════════════════════════╝

def make_summary_tsv(df, refs):
    """Build detailed summary output using evidence-oriented display fields."""

    def normalize_scores(dataf, sample_col='sample_id', score_col='predicted_score'):
        dataf = dataf.copy()
        for sample in dataf[sample_col].unique():
            mask = dataf[sample_col] == sample
            min_s = dataf.loc[mask, score_col].min()
            max_s = dataf.loc[mask, score_col].max()
            dataf.loc[mask, 'normalized_score'] = (
                (dataf.loc[mask, score_col] - min_s) / (max_s - min_s)
                if max_s > min_s else 1.0
            )
        return dataf

    def build_summary_lookups(refs):
        drug_df = pd.read_csv(refs['druggable_file'], sep='	')
        panel_df = pd.read_csv(refs['fusion_panels'], sep='	')
        cancermine_df = pd.read_csv(refs['cancermine_file'], sep='	')

        cm_pivot = cancermine_df.groupby(['gene_normalized', 'role'])['citation_count'].sum().unstack(fill_value=0)
        cm_pivot.columns = [c.upper() for c in cm_pivot.columns]
        cm_pivot = cm_pivot.rename(columns={'TUMOR_SUPPRESSOR': 'TSG'})
        cm_pivot.index = cm_pivot.index.str.upper()

        panel_cols = ['chop', 'stjude', 'oncokids', 'archer', 'caris',
                      'msk_impact', 'msk_heme', 'foundation', 'foundation_heme',
                      'vogelstein', 'oncokb']
        gene_to_panels = {}
        for row in panel_df.itertuples(index=False):
            gene = str(getattr(row, 'gene')).upper()
            panels = [col for col in panel_cols if getattr(row, col, 0) == 1]
            if panels:
                gene_to_panels[gene] = panels

        fusion_to_drugs = {}
        gene_to_drug_cancers = {}
        for row in drug_df.itertuples(index=False):
            g1 = str(getattr(row, 'gene1')).upper()
            g2 = str(getattr(row, 'gene2')).upper()
            gene_t = str(getattr(row, 'gene_targets')).upper()
            cancers = str(getattr(row, 'all_actionable_cancers')).upper()
            entry = (gene_t, cancers)
            fusion_to_drugs.setdefault((g1, g2), []).append(entry)
            fusion_to_drugs.setdefault((g2, g1), []).append(entry)
            gene_to_drug_cancers.setdefault(gene_t, set()).add(cancers)

        return cm_pivot, gene_to_panels, fusion_to_drugs, gene_to_drug_cancers

    def get_gene_roles(gene, is_kinase, is_rtk, is_tf, is_onco, is_tsg, cm_pivot):
        gene_upper = str(gene).upper()
        roles = []

        if is_rtk:
            roles.append('RTK')
        elif is_kinase:
            roles.append('kinase')
        if is_tf:
            roles.append('TF')

        onco_cit = cm_pivot.loc[gene_upper, 'ONCOGENE'] if gene_upper in cm_pivot.index and 'ONCOGENE' in cm_pivot.columns else 0
        tsg_cit = cm_pivot.loc[gene_upper, 'TSG'] if gene_upper in cm_pivot.index and 'TSG' in cm_pivot.columns else 0
        total = onco_cit + tsg_cit

        if is_kinase or is_rtk:
            if is_onco or is_tsg:
                roles.append('oncogene')
            if is_tsg and total > 50 and tsg_cit / total > 0.90:
                roles.append('TSG (strong evidence)')
        elif is_onco and is_tsg:
            if total == 0:
                roles.append('oncogene')
            elif onco_cit >= tsg_cit:
                roles.append('oncogene')
                if tsg_cit > 0:
                    roles.append('TSG (conflicting)')
            else:
                roles.append('TSG')
                if onco_cit > 0:
                    roles.append('oncogene (conflicting)')
        elif is_onco:
            roles.append('oncogene')
        elif is_tsg:
            if total > 0 and tsg_cit / total > 0.75:
                roles.append('TSG')
            else:
                roles.append('TSG (low evidence)')
        return roles

    df = normalize_scores(df)
    cm_pivot, gene_to_panels, fusion_to_drugs, gene_to_drug_cancers = build_summary_lookups(refs)

    key_cols = ['sample_id', 'gene1', 'gene2']
    work = df.sort_values(['sample_id', 'gene1', 'gene2', 'normalized_score'], ascending=[True, True, True, False]).copy()

    # Precompute all breakpoints text + best breakpoint row once per fusion key.
    has_bkp = all(c in work.columns for c in ['gene1_chr', 'gene1_bkp', 'gene2_chr', 'gene2_bkp'])
    if has_bkp:
        bkp_cols = ['sample_id', 'gene1', 'gene2', 'gene1_chr', 'gene1_bkp', 'gene2_chr', 'gene2_bkp', 'normalized_score']
        bkp_df = work[bkp_cols].drop_duplicates()
    else:
        bkp_cols = ['sample_id', 'gene1', 'gene2', 'normalized_score']
        bkp_df = work[bkp_cols].drop_duplicates()


        
    def _bkp_text(g):
        if has_bkp:
            return '\n'.join(
                f"{g.name[1]};{r.gene1_chr}:{int(r.gene1_bkp)},{g.name[2]};{r.gene2_chr}:{int(r.gene2_bkp)};{r.normalized_score:.3f}"
                for r in g.itertuples(index=False)
            )
        else:
            return '\n'.join(
                f"{g.name[1]}::{g.name[2]};{r.normalized_score:.3f}"
                for r in g.itertuples(index=False)
            )        

    bkp_text = (
        bkp_df.groupby(key_cols, sort=False)
        .apply(_bkp_text, include_groups=False)
        .rename('all_breakpoints')
        .reset_index()
    )

    best_bkp = bkp_df.drop_duplicates(subset=key_cols, keep='first').copy()
    
    if has_bkp:
        best_bkp['gene1_breakpoint'] = best_bkp.apply(lambda r: f"{r['gene1_chr']}:{int(r['gene1_bkp']):,}", axis=1)
        best_bkp['gene2_breakpoint'] = best_bkp.apply(lambda r: f"{r['gene2_chr']}:{int(r['gene2_bkp']):,}", axis=1)
    else:
        best_bkp['gene1_breakpoint'] = '—'
        best_bkp['gene2_breakpoint'] = '—'    
    
    best_bkp = best_bkp[key_cols + ['gene1_breakpoint', 'gene2_breakpoint']]

    sample_data = work.drop_duplicates(subset=key_cols, keep='first').copy()
    sample_data['Rank'] = sample_data.groupby('sample_id')['normalized_score'].rank(ascending=False, method='first').astype(int)

    sample_data = sample_data.merge(best_bkp, on=key_cols, how='left').merge(bkp_text, on=key_cols, how='left')

    rows = []
    for row in sample_data.itertuples(index=False):
        db_evidence = []
        if getattr(row, 'in_civic_exact', 0):
            db_evidence.append('CIViC exact match')
        if getattr(row, 'in_mitelman_exact', 0):
            db_evidence.append('Mitelman DB match')
        if getattr(row, 'in_cosmic_exact', 0):
            db_evidence.append('COSMIC exact match')
        if getattr(row, 'in_cosmic_same_cancer', 0):
            db_evidence.append('COSMIC same cancer')
        if getattr(row, 'cosmic_gene_involved', 0) and not getattr(row, 'in_cosmic_exact', 0):
            db_evidence.append('COSMIC gene involved')
        db_str = '\n'.join(db_evidence) if db_evidence else '—'

        drug_evidence = []
        if getattr(row, 'civic_has_therapeutic', 0):
            drug_evidence.append('Approved therapy (CIViC)')

        ct = str(getattr(row, 'cancer_type', ''))
        ct_upper = ct.upper()
        fkey = (str(getattr(row, 'gene1')).upper(), str(getattr(row, 'gene2')).upper())
        if fkey in fusion_to_drugs:
            cancers = {c for _, c in fusion_to_drugs[fkey]}
            if ct_upper in cancers:
                drug_evidence.append(f"{row.gene1}::{row.gene2} — druggable in {ct}")
            else:
                drug_evidence.append(f"{row.gene1}::{row.gene2} — druggable (other cancers)")
        else:
            for gene in [str(row.gene1).upper(), str(row.gene2).upper()]:
                if gene in gene_to_drug_cancers:
                    cancers = gene_to_drug_cancers[gene]
                    if ct_upper in cancers:
                        drug_evidence.append(f"{gene} — druggable in {ct}")
                    else:
                        drug_evidence.append(f"{gene} — druggable (other cancers)")
        drug_str = '\n'.join(drug_evidence) if drug_evidence else '—'

        biology = []
        for gene, is_rtk, is_kinase, is_tf, is_onco, is_tsg in [
            (row.gene1, getattr(row, 'gene1_is_rtk', 0), getattr(row, 'gene1_is_kinase', 0),
             getattr(row, 'gene1_is_tf', 0), getattr(row, 'gene1_oncogene', 0), getattr(row, 'gene1_tsg', 0)),
            (row.gene2, getattr(row, 'gene2_is_rtk', 0), getattr(row, 'gene2_is_kinase', 0),
             getattr(row, 'gene2_is_tf', 0), getattr(row, 'gene2_oncogene', 0), getattr(row, 'gene2_tsg', 0)),
        ]:
            roles = get_gene_roles(gene, is_kinase, is_rtk, is_tf, is_onco, is_tsg, cm_pivot)
            if roles:
                biology.append(f"{gene}: {', '.join(roles)}")
        bio_str = '\n'.join(biology) if biology else '—'

        panel_genes = [g for g in [str(row.gene1).upper(), str(row.gene2).upper()] if g in gene_to_panels]
        panel_str = ', '.join(panel_genes) if panel_genes else '—'

        rows.append({
            'Sample': row.sample_id,
            'Cancer': ct,
            'Rank': int(row.Rank),
            'Fusion': row.fusion_name,
            'Score (0-1)': round(float(row.normalized_score), 3),
            'Gene1 Breakpoint': row.gene1_breakpoint if pd.notna(row.gene1_breakpoint) else '—',
            'Gene2 Breakpoint': row.gene2_breakpoint if pd.notna(row.gene2_breakpoint) else '—',
            'Database Evidence': db_str,
            'Drug Evidence': drug_str,
            'Gene Biology': bio_str,
            'Panel Evidence': panel_str,
            'All Breakpoints (Norm. Score)': row.all_breakpoints if pd.notna(row.all_breakpoints) else '—',
        })

    result = pd.DataFrame(rows)
    display_cols = ['Sample', 'Cancer', 'Rank', 'Fusion', 'Score (0-1)',
                    'Gene1 Breakpoint', 'Gene2 Breakpoint', 'Database Evidence',
                    'Drug Evidence', 'Gene Biology', 'Panel Evidence',
                    'All Breakpoints (Norm. Score)']
    if len(result):
        result = result.sort_values(['Sample', 'Rank']).reset_index(drop=True)
        return result[display_cols]
    return pd.DataFrame(columns=display_cols)


def write_summary_html(summary_df, output_path, sample_prefix):
    """Write an easy-to-read HTML report for the ranked fusion summary."""

    def fmt(value):
        if pd.isna(value):
            return '—'
        return escape(str(value)).replace('\n', '<br>')

    def row_class(row):
        try:
            rank = int(row['Rank'])
        except (TypeError, ValueError):
            rank = None
        return ' class="top-rank"' if rank is not None and rank <= 5 else ''

    columns = list(summary_df.columns)
    generated_label = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
    row_count = len(summary_df)
    sample_count = summary_df['Sample'].nunique() if 'Sample' in summary_df.columns and row_count else 0

    header_html = ''.join(f'<th>{escape(str(col))}</th>' for col in columns)
    if row_count:
        body_rows = []
        for _, row in summary_df.iterrows():
            cells = ''.join(f'<td>{fmt(row[col])}</td>' for col in columns)
            body_rows.append(f'<tr{row_class(row)}>{cells}</tr>')
        body_html = '\n'.join(body_rows)
    else:
        colspan = max(len(columns), 1)
        body_html = f'<tr><td class="empty" colspan="{colspan}">No fusion candidates were reported.</td></tr>'

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>FusionRank Summary - {escape(str(sample_prefix))}</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f8fb;
      --panel: #ffffff;
      --text: #172033;
      --muted: #5f6b7a;
      --border: #d9e2ec;
      --header: #12355b;
      --header-text: #ffffff;
      --accent: #e8f3ff;
      --top: #f0fdf4;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      line-height: 1.45;
    }}
    main {{
      max-width: 1400px;
      margin: 0 auto;
      padding: 32px 24px;
    }}
    .report-header {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 14px;
      box-shadow: 0 8px 24px rgba(15, 23, 42, 0.06);
      margin-bottom: 24px;
      padding: 24px;
    }}
    h1 {{
      margin: 0 0 8px;
      color: var(--header);
      font-size: 1.9rem;
    }}
    .meta {{
      color: var(--muted);
      display: flex;
      flex-wrap: wrap;
      gap: 16px;
      margin-top: 12px;
    }}
    .meta strong {{ color: var(--text); }}
    .table-wrap {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 14px;
      box-shadow: 0 8px 24px rgba(15, 23, 42, 0.06);
      overflow: auto;
    }}
    table {{
      border-collapse: separate;
      border-spacing: 0;
      min-width: 1200px;
      width: 100%;
    }}
    th, td {{
      border-bottom: 1px solid var(--border);
      padding: 12px 14px;
      text-align: left;
      vertical-align: top;
      white-space: normal;
    }}
    th {{
      background: var(--header);
      color: var(--header-text);
      font-size: 0.78rem;
      letter-spacing: 0.03em;
      position: sticky;
      text-transform: uppercase;
      top: 0;
      z-index: 1;
    }}
    tbody tr:nth-child(even) {{ background: #fbfdff; }}
    tbody tr:hover {{ background: var(--accent); }}
    tbody tr.top-rank {{ background: var(--top); }}
    td:nth-child(3), td:nth-child(5) {{
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }}
    td:nth-child(4) {{ font-weight: 700; }}
    .empty {{
      color: var(--muted);
      text-align: center;
    }}
    .hint {{
      color: var(--muted);
      margin: 16px 0 0;
    }}
  </style>
</head>
<body>
  <main>
    <section class="report-header">
      <h1>FusionRank Summary</h1>
      <p class="hint">Ranked fusion candidates with evidence fields formatted for browser review. Rows ranked 1–5 are highlighted.</p>
      <div class="meta">
        <span><strong>Sample prefix:</strong> {escape(str(sample_prefix))}</span>
        <span><strong>Samples:</strong> {sample_count:,}</span>
        <span><strong>Fusions:</strong> {row_count:,}</span>
        <span><strong>Generated:</strong> {generated_label}</span>
      </div>
    </section>
    <section class="table-wrap" aria-label="FusionRank summary table">
      <table>
        <thead><tr>{header_html}</tr></thead>
        <tbody>
{body_html}
        </tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    with open(output_path, 'w', encoding='utf-8') as handle:
        handle.write(html_text)



# ANNOTATE ONLY ---- BEFORE
# df = run_full_pipeline(input_file=<input_file>, output_dir=<output_dir>, project_name=<project_name>, refs=<refs>);

# ╔══════════════════════════════════════════════════════════════╗
# ║  RUN:   ANNOTATION - PREDICTED RANKING - SUMMARY OUTPUTS     ║
# ╚══════════════════════════════════════════════════════════════╝

