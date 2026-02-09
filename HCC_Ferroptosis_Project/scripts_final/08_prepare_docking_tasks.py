#!/usr/bin/env python3

"""
Prepare CB-Dock2 task list from Ezhu ligands with a reproducible filter.

Selection rules:
1) Ligand name contains letters (avoid pure CAS/number labels).
2) SMILES is present.
3) Deduplicate by SMILES (first occurrence kept).
4) Select first N ligands by original order to avoid manual bias.
"""

import argparse
import csv
import os
import re


def load_ligands(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    return reader.fieldnames, rows


def filter_ligands(rows, max_n):
    seen_smiles = set()
    selected = []
    for row in rows:
        name = (row.get("Ligand") or "").strip()
        smiles = (row.get("SMILES") or "").strip()
        if not name or not smiles:
            continue
        if not re.search(r"[A-Za-z]", name):
            continue
        if smiles in seen_smiles:
            continue
        seen_smiles.add(smiles)
        selected.append(row)
        if len(selected) >= max_n:
            break
    return selected


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Prepare docking task list for CB-Dock2.")
    parser.add_argument("--ligands", default="data/references/docking/ligands_smiles_ezhu.csv")
    parser.add_argument("--targets", default="AKT1,PIK3CA,MTOR")
    parser.add_argument("--max-ligands", type=int, default=15)
    parser.add_argument("--out-ligands", default="data/references/docking/ligands_smiles_ezhu_top15.csv")
    parser.add_argument("--out-tasks", default="results/docking_tasks.csv")
    args = parser.parse_args()

    fieldnames, rows = load_ligands(args.ligands)
    selected = filter_ligands(rows, args.max_ligands)
    if not selected:
        raise SystemExit("No ligands selected. Check ligand file and filters.")

    os.makedirs(os.path.dirname(args.out_ligands), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_tasks), exist_ok=True)

    write_csv(args.out_ligands, fieldnames, selected)

    targets = [t.strip().upper() for t in args.targets.split(",") if t.strip()]
    task_rows = []
    for t in targets:
        for row in selected:
            task_rows.append(
                {
                    "Target": t,
                    "Ligand": row.get("Ligand", "").strip(),
                    "SMILES": row.get("SMILES", "").strip(),
                }
            )

    write_csv(args.out_tasks, ["Target", "Ligand", "SMILES"], task_rows)

    print(f"Selected ligands: {len(selected)}")
    print(f"Targets: {', '.join(targets)}")
    print(f"Tasks: {len(task_rows)}")
    print(f"Output ligands: {args.out_ligands}")
    print(f"Output tasks: {args.out_tasks}")


if __name__ == "__main__":
    main()
