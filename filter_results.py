import os
import pandas as pd
import re
import streamlit as st

def clean_sheet_name(name):
    base_name = os.path.basename(name)
    clean_name = re.sub(r'[\\/*?:\[\]]', '_', base_name.split('.')[0])
    return clean_name[:31]

def process_rmats_output(output_dir, pval_cutoff):
    rmats_out = os.path.join(output_dir, "output_combined")
    files_to_check = ['A3SS.MATS.JCEC.txt', 'A5SS.MATS.JCEC.txt', 'MXE.MATS.JCEC.txt', 'RI.MATS.JCEC.txt', 'SE.MATS.JCEC.txt']
    
    existing_files = [os.path.join(rmats_out, f) for f in files_to_check if os.path.exists(os.path.join(rmats_out, f))]
    
    if not existing_files:
        return None, 0
        
    all_results = {}
    for file_path in existing_files:
        try:
            df = pd.read_csv(file_path, sep='\t')
            
            if 'longExonStart_0base' in df.columns:
                df['UCSC_Coordinate'] = df['chr'] + ":" + (df['longExonStart_0base'] + 1).astype(str) + "-" + df['longExonEnd'].astype(str)
            elif 'exonStart_0base' in df.columns:
                df['UCSC_Coordinate'] = df['chr'] + ":" + (df['exonStart_0base'] + 1).astype(str) + "-" + df['exonEnd'].astype(str)
            elif 'riExonStart_0base' in df.columns:
                df['UCSC_Coordinate'] = df['chr'] + ":" + (df['riExonStart_0base'] + 1).astype(str) + "-" + df['riExonEnd'].astype(str)
            elif '1stExonStart_0base' in df.columns:
                df['UCSC_Exon1'] = df['chr'] + ":" + (df['1stExonStart_0base'] + 1).astype(str) + "-" + df['1stExonEnd'].astype(str)
                df['UCSC_Exon2'] = df['chr'] + ":" + (df['2ndExonStart_0base'] + 1).astype(str) + "-" + df['2ndExonEnd'].astype(str)
            
            def get_ijc_avg(val):
                if pd.isna(val) or val == "NA": return 0
                nums = [float(i) for i in str(val).split(',')]
                return sum(nums) / len(nums)

            df['avg_IJC1'] = df['IJC_SAMPLE_1'].apply(get_ijc_avg)
            df['avg_IJC2'] = df['IJC_SAMPLE_2'].apply(get_ijc_avg)
            
            # Apply dynamic PValue/FDR filter from UI
            coverage_mask = ((df['avg_IJC1'] + df['avg_IJC2']) / 2) >= 20
            significance_mask = (df['PValue'] < pval_cutoff) & (df['FDR'] < pval_cutoff)
            
            filtered_df = df[coverage_mask & significance_mask].copy()
            
            cols_to_keep = ['ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'PValue', 'FDR', 'IncLevelDifference']
            coord_cols = [c for c in filtered_df.columns if 'UCSC' in c]
            final_df = filtered_df[cols_to_keep + coord_cols]
            
            all_results[clean_sheet_name(file_path)] = final_df
            
        except Exception as e:
            st.error(f"Error processing {file_path}: {e}")

    if all_results:
        filter_out_dir = os.path.join(output_dir, 'filtered_results')
        os.makedirs(filter_out_dir, exist_ok=True)
        output_xlsx = os.path.join(filter_out_dir, 'Filtered_Splicing_Results.xlsx')
        
        with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
            for sheet_name, data in all_results.items():
                data.to_excel(writer, sheet_name=sheet_name, index=False)
                
        total_events = sum(len(v) for v in all_results.values())
        return output_xlsx, total_events
    
    return None, 0