#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 21:15:47 2025

@author: chenmenghui
"""

import pandas as pd
import glob, os, re

INPUT_DIR   = "/Users/chenmenghui/Desktop/shinyProject/datasets/pval_results/"
OUTPUT_XLSX = os.path.join(INPUT_DIR, "age_sig_genes_by_tissue.xlsx")

PVAL_COL    = "BH_adjusted_age"
ALPHA       = 0.05              
DROP_NA     = True                
WRITE_SUMMARY = True           

def find_pval_col(cols, target="BH_adjusted_age"):
    for c in cols:
        if c.strip().lower() == target.lower():
            return c
    cand = [c for c in cols if ("bh" in c.lower() and "adjust" in c.lower() and "age" in c.lower())]
    return cand[0] if cand else None

def clean_sheet_name(name, used):
    for suf in ["_pvalue_results", "-pvalue-results"]:
        if name.endswith(suf):
            name = name[: -len(suf)]
    name = re.sub(r'[\[\]\*\?\\/:\']', '_', name).strip()
    if len(name) > 31:
        name = name[:31]
    base = name
    i = 1
    while name in used:
        suffix = f"_{i}"
        allowed = 31 - len(suffix)
        name = (base[:allowed] + suffix) if len(base) + len(suffix) > 31 else base + suffix
        i += 1
    used.add(name)
    return name

def main():
    csv_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.csv")))
    if not csv_files:
        print("No CSV files were found. Please check INPUT_DIR.")
        return

    sheet_names_used = set()
    summary_rows = []
    sheets_written = 0

    with pd.ExcelWriter(OUTPUT_XLSX, engine="xlsxwriter") as writer:
        for f in csv_files:
            try:
                df = pd.read_csv(f)
            except Exception as e:
                summary_rows.append({"File": os.path.basename(f), "Status": f"读入失败: {e}", 
                                     "Rows_total": None, "Rows_kept": None})
                continue

            df.columns = [c.strip() for c in df.columns]
            pcol = find_pval_col(df.columns, PVAL_COL)
            if pcol is None:
                summary_rows.append({"File": os.path.basename(f), "Status": "缺少 BH_adjusted_age 列", 
                                     "Rows_total": int(df.shape[0]), "Rows_kept": 0})
                continue

            # 转数值并筛选
            df["_bh_num"] = pd.to_numeric(df[pcol], errors="coerce")
            filtered = df[df["_bh_num"] < ALPHA].copy()
            if DROP_NA:
                filtered = filtered.dropna(axis=0, how="any")
            filtered.drop(columns=["_bh_num"], inplace=True)

            kept = int(filtered.shape[0])
            total = int(df.shape[0])

            # 只写有数据的 sheet
            if kept > 0:
                base = os.path.splitext(os.path.basename(f))[0]
                sheet = clean_sheet_name(base, sheet_names_used)
                filtered.to_excel(writer, index=False, sheet_name=sheet)
                sheets_written += 1
                summary_rows.append({"File": os.path.basename(f), "Status": "OK", 
                                     "Rows_total": total, "Rows_kept": kept, "Sheet": sheet})
            else:
                summary_rows.append({"File": os.path.basename(f), "Status": "No age-related significant genes (no table established)", 
                                     "Rows_total": total, "Rows_kept": kept})

        if WRITE_SUMMARY:
            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_excel(writer, index=False, sheet_name="Summary")

    print(f"Completed: Write {sheets_written} sheets containing data → {OUTPUT_XLSX}")

if __name__ == "__main__":
    main()