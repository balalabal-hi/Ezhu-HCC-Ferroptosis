#!/usr/bin/env python3

"""
Clean CB-Dock2 results and map ligand names back to canonical labels.

Input:
  data/references/docking/docking_results.csv
  data/references/docking/ligands_sdf/ligands_sdf_map.csv

Outputs:
  results/cbdock2_docking_results.csv
  results/cbdock2_top_pairs.csv
"""

import csv
import os


def load_map(path):
    mapping = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sdf_file = (row.get("SDF_File") or "").strip()
            ligand = (row.get("Ligand") or "").strip()
            if sdf_file and ligand:
                base = os.path.splitext(sdf_file)[0]
                mapping[base] = ligand
    return mapping


def to_float(val):
    try:
        return float(val)
    except Exception:
        return None


def main():
    in_results = "data/references/docking/docking_results.csv"
    map_path = "data/references/docking/ligands_sdf/ligands_sdf_map.csv"
    out_results = "results/cbdock2_docking_results.csv"
    out_top = "results/cbdock2_top_pairs.csv"

    if not os.path.exists(in_results):
        raise SystemExit(f"Missing {in_results}")
    if not os.path.exists(map_path):
        raise SystemExit(f"Missing {map_path}")

    ligand_map = load_map(map_path)

    cleaned = []
    with open(in_results, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ligand_raw = (row.get("Ligand") or "").strip()
            ligand = ligand_map.get(ligand_raw, ligand_raw)
            contact = (row.get("Contact_Residues") or "").strip()
            if contact.lower() == "view":
                contact = ""

            docking_size_x = to_float(row.get("Docking_Size_X"))
            docking_size_y = to_float(row.get("Docking_Size_Y"))
            docking_size_z = to_float(row.get("Docking_Size_Z"))
            valid_size = all(
                v is not None and v > 0 for v in [docking_size_x, docking_size_y, docking_size_z]
            )
            if not valid_size:
                continue

            cleaned.append(
                {
                    "Protein": (row.get("Protein") or "").strip(),
                    "Ligand": ligand,
                    "Vina_Score": to_float(row.get("Vina_Score")),
                    "Cavity_Volume": to_float(row.get("Cavity_Volume")),
                    "Center_X": to_float(row.get("Center_X")),
                    "Center_Y": to_float(row.get("Center_Y")),
                    "Center_Z": to_float(row.get("Center_Z")),
                    "Docking_Size_X": docking_size_x,
                    "Docking_Size_Y": docking_size_y,
                    "Docking_Size_Z": docking_size_z,
                    "Contact_Residues": contact,
                }
            )

    os.makedirs(os.path.dirname(out_results), exist_ok=True)

    with open(out_results, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Protein",
                "Ligand",
                "Vina_Score",
                "Cavity_Volume",
                "Center_X",
                "Center_Y",
                "Center_Z",
                "Docking_Size_X",
                "Docking_Size_Y",
                "Docking_Size_Z",
                "Contact_Residues",
            ],
        )
        writer.writeheader()
        for row in cleaned:
            writer.writerow(row)

    # Top pairs by protein
    top_rows = []
    by_protein = {}
    for row in cleaned:
        prot = row["Protein"]
        if not prot:
            continue
        by_protein.setdefault(prot, []).append(row)

    for prot, rows in by_protein.items():
        rows = [r for r in rows if r["Vina_Score"] is not None]
        if not rows:
            continue
        best = min(rows, key=lambda r: r["Vina_Score"])
        best = best.copy()
        best["Rank"] = 1
        top_rows.append(best)

    with open(out_top, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Protein",
                "Ligand",
                "Vina_Score",
                "Cavity_Volume",
                "Center_X",
                "Center_Y",
                "Center_Z",
                "Docking_Size_X",
                "Docking_Size_Y",
                "Docking_Size_Z",
                "Contact_Residues",
                "Rank",
            ],
        )
        writer.writeheader()
        for row in top_rows:
            writer.writerow(row)

    print(f"Wrote {out_results}")
    print(f"Wrote {out_top}")


if __name__ == "__main__":
    main()
