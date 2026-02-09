#!/usr/bin/env python3

"""
Run Boltz-2 NIM affinity predictions for docking pairs.

Inputs:
  - results/cbdock2_docking_results.csv (default, all pairs)
  - results/cbdock2_top_pairs.csv (optional, top-only)
  - data/references/docking/ligands_smiles_ezhu_top15.csv (or full list)

Outputs:
  - results/boltz2_affinity_results.csv
  - results/boltz2_responses/*.json (raw responses)

Auth:
  - Set NVIDIA_API_KEY in environment (do not hardcode; do not include keys in shared packages).

Endpoint:
  - Use BOLTZ2_ENDPOINTS to override, comma-separated.
  - Default attempts:
      https://health.api.nvidia.com/v1/biology/mit/boltz2/predict
      https://api.nvidia.com/v1/biology/mit/boltz2/predict
      https://integrate.api.nvidia.com/v1/biology/mit/boltz2/predict
"""

import csv
import json
import os
import re
import time
from urllib.request import urlopen

import argparse
try:
    import requests  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    requests = None


UNIPROT_IDS = {
    "AKT1": "P31749",
    "PIK3CA": "P42336",
    "MTOR": "P42345",
    "MAPK8": "P45983",
    "SLC27A5": "Q9Y2P5",
    "ENO1": "P06733",
}

STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{task_id}"
MTOR_DOMAIN_RANGE = os.getenv("MTOR_DOMAIN_RANGE", "2156-2469")


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_")


def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    with urlopen(url) as resp:
        data = resp.read().decode("utf-8")
    seq = "".join(line.strip() for line in data.splitlines() if not line.startswith(">"))
    if not seq:
        raise RuntimeError(f"Empty sequence for {uniprot_id}")
    return seq


def extract_sequence_from_pdb(pdb_path, chain_id="A"):
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54:
                continue
            if line[21].strip() != chain_id:
                continue
            res_name = line[17:20].strip()
            try:
                res_id = int(line[22:26].strip())
            except ValueError:
                continue
            residues[res_id] = res_name
    if not residues:
        raise RuntimeError(f"No residues parsed from {pdb_path} chain {chain_id}")
    aa_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    ordered = [aa_map.get(residues[r], "X") for r in sorted(residues)]
    return "".join(ordered), sorted(residues)


def get_mtor_domain_sequence():
    pdb_path = os.path.join("data", "references", "docking", "targets", "MTOR.pdb")
    if not os.path.exists(pdb_path):
        raise RuntimeError("MTOR.pdb not found for domain extraction.")
    full_seq, res_ids = extract_sequence_from_pdb(pdb_path, chain_id="A")
    start, end = MTOR_DOMAIN_RANGE.split("-")
    start = int(start)
    end = int(end)
    if start not in res_ids or end not in res_ids:
        raise RuntimeError("MTOR domain range not present in MTOR.pdb residues.")
    # Map residue numbers to sequence positions
    idx_map = {res_id: i for i, res_id in enumerate(res_ids)}
    s_idx = idx_map[start]
    e_idx = idx_map[end]
    return full_seq[s_idx : e_idx + 1]


def load_smiles_map(paths):
    mapping = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                ligand = (row.get("Ligand") or "").strip()
                smiles = (row.get("SMILES") or "").strip()
                if ligand and smiles and ligand not in mapping:
                    mapping[ligand] = smiles
    return mapping


def normalize_ligand(name):
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def build_normalized_map(smiles_map):
    normalized = {}
    for ligand, smiles in smiles_map.items():
        key = normalize_ligand(ligand)
        if key and key not in normalized:
            normalized[key] = smiles
    return normalized


def get_endpoints():
    env = os.getenv("BOLTZ2_ENDPOINTS", "").strip()
    if env:
        return [e.strip() for e in env.split(",") if e.strip()]
    return [
        "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict",
        "https://api.nvidia.com/v1/biology/mit/boltz2/predict",
        "https://integrate.api.nvidia.com/v1/biology/mit/boltz2/predict",
    ]


def call_boltz2(api_key, endpoint, sequence, smiles, contacts=None, auth_mode="bearer"):
    if requests is None:
        raise RuntimeError("Python package 'requests' is required for Boltz-2 API calls. Install it or use --collect-only.")
    headers = {
        "Content-Type": "application/json",
        "accept": "application/json",
    }
    if auth_mode == "bearer":
        headers["Authorization"] = f"Bearer {api_key}"
    elif auth_mode == "nvidia":
        headers["NVIDIA-API-Key"] = api_key
    else:
        raise ValueError("auth_mode must be 'bearer' or 'nvidia'")
    headers["NVCF-POLL-SECONDS"] = os.getenv("NVCF_POLL_SECONDS", "300")
    payload = {
        "polymers": [
            {
                "id": "A",
                "molecule_type": "protein",
                "sequence": sequence,
            }
        ],
        "ligands": [
            {
                "id": "L1",
                "smiles": smiles,
                "predict_affinity": True,
            }
        ],
        "sampling_steps": int(os.getenv("BOLTZ2_SAMPLING_STEPS", "50")),
        "diffusion_samples": int(os.getenv("BOLTZ2_DIFFUSION_SAMPLES", "1")),
        "sampling_steps_affinity": int(os.getenv("BOLTZ2_SAMPLING_STEPS_AFFINITY", "50")),
        "diffusion_samples_affinity": int(os.getenv("BOLTZ2_DIFFUSION_SAMPLES_AFFINITY", "1")),
        "output_format": "mmcif",
    }
    if contacts:
        payload["constraints"] = [
            {
                "constraint_type": "pocket",
                "binder": "L1",
                "contacts": contacts,
            }
        ]
    resp = requests.post(
        endpoint,
        headers=headers,
        data=json.dumps(payload),
        timeout=float(os.getenv("BOLTZ2_TIMEOUT", "400")),
    )
    if resp.status_code != 202:
        return resp

    task_id = resp.headers.get("nvcf-reqid")
    if not task_id:
        return resp

    poll_headers = {
        "Authorization": headers.get("Authorization", ""),
        "NVIDIA-API-Key": headers.get("NVIDIA-API-Key", ""),
        "accept": "application/json",
    }
    poll_headers = {k: v for k, v in poll_headers.items() if v}
    poll_interval = float(os.getenv("NVCF_POLL_INTERVAL", "10"))
    max_wait = float(os.getenv("NVCF_POLL_TIMEOUT", "900"))
    start = time.time()
    while time.time() - start < max_wait:
        status = requests.get(
            STATUS_URL.format(task_id=task_id),
            headers=poll_headers,
            timeout=float(os.getenv("BOLTZ2_TIMEOUT", "400")),
        )
        if status.status_code == 200:
            return status
        if status.status_code in {202, 404}:
            time.sleep(poll_interval)
            continue
        if status.status_code in {400, 401, 403, 422, 500}:
            return status
        time.sleep(poll_interval)
    return resp


def parse_affinity(result):
    affinity_val = None
    affinity_prob = None
    affinity_pic50 = None
    confidence = None

    if isinstance(result, dict):
        affinity_val = result.get("affinity_pred_value") or result.get("affinity_value")
        affinity_prob = result.get("affinity_probability_binary") or result.get("affinity_probability")
        conf = result.get("confidence_scores")
        if isinstance(conf, list) and conf:
            confidence = conf[0]
        affinities = result.get("affinities")
        if affinity_val is None and isinstance(affinities, dict) and affinities:
            affinity_entry = None
            if "L1" in affinities and isinstance(affinities["L1"], dict):
                affinity_entry = affinities["L1"]
            else:
                affinity_entry = next(iter(affinities.values()))
            if isinstance(affinity_entry, dict):
                affinity_val = affinity_entry.get("affinity_pred_value") or affinity_entry.get("affinity_value")
                affinity_prob = affinity_entry.get("affinity_probability_binary") or affinity_entry.get("affinity_probability")
                affinity_pic50 = affinity_entry.get("affinity_pic50")
        elif affinity_val is None and isinstance(affinities, list) and affinities:
            first = affinities[0]
            if isinstance(first, dict):
                affinity_val = first.get("affinity_pred_value") or first.get("affinity_value")
                affinity_prob = first.get("affinity_probability_binary") or first.get("affinity_probability")
                affinity_pic50 = first.get("affinity_pic50")
        elif isinstance(affinities, dict) and "L1" in affinities and isinstance(affinities["L1"], dict):
            affinity_pic50 = affinities["L1"].get("affinity_pic50")
    if isinstance(affinity_val, list) and affinity_val:
        affinity_val = affinity_val[0]
    if isinstance(affinity_prob, list) and affinity_prob:
        affinity_prob = affinity_prob[0]
    if isinstance(affinity_pic50, list) and affinity_pic50:
        affinity_pic50 = affinity_pic50[0]
    return affinity_val, affinity_pic50, affinity_prob, confidence


def parse_contacts(text, max_contacts=20):
    if not text:
        return []
    tokens = [t.strip() for t in text.split(",") if t.strip()]
    contacts = []
    for tok in tokens:
        parts = tok.split(":")
        if len(parts) < 2:
            continue
        # Accept formats like GLU:17:A or GLU:17
        residue_index = None
        chain_id = "A"
        if len(parts) >= 3:
            try:
                residue_index = int(parts[1])
                chain_id = parts[2] or "A"
            except Exception:
                residue_index = None
        else:
            try:
                residue_index = int(parts[1])
            except Exception:
                residue_index = None
        if residue_index is None:
            continue
        contacts.append({"id": chain_id, "residue_index": residue_index})
        if len(contacts) >= max_contacts:
            break
    return contacts


def main():
    parser = argparse.ArgumentParser(description="Boltz-2 affinity scoring for docking pairs.")
    parser.add_argument("--pairs", default="results/cbdock2_docking_results.csv")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no-resume", dest="resume", action="store_false")
    parser.add_argument("--max", type=int, default=0, help="Limit number of pairs (0 = all).")
    parser.add_argument("--only-protein", default="", help="Process only this protein (exact match).")
    parser.add_argument("--only-ligand", default="", help="Process only this ligand (exact match).")
    parser.add_argument("--use-pocket", action="store_true", default=True)
    parser.add_argument("--no-pocket", dest="use_pocket", action="store_false")
    parser.add_argument("--sleep", type=float, default=1.0)
    parser.add_argument("--test", action="store_true", default=False, help="Test a single request and exit.")
    parser.add_argument("--collect-only", action="store_true", default=False, help="Only aggregate existing responses.")
    args = parser.parse_args()

    pairs_path = args.pairs
    if not os.path.exists(pairs_path):
        raise SystemExit(f"Missing {pairs_path}")

    smiles_map = load_smiles_map(
        [
            "data/references/docking/ligands_smiles_ezhu_top15.csv",
            "data/references/docking/ligands_smiles_ezhu.csv",
        ]
    )
    if not smiles_map:
        raise SystemExit("Missing ligand SMILES mapping.")
    normalized_smiles_map = build_normalized_map(smiles_map)

    os.makedirs("results/boltz2_responses", exist_ok=True)
    os.makedirs("results/boltz2_failures", exist_ok=True)

    endpoints = get_endpoints()

    with open(pairs_path, newline="") as f:
        pairs = list(csv.DictReader(f))

    if args.max and args.max > 0:
        pairs = pairs[: args.max]

    if args.collect_only:
        results = []
        for row in pairs:
            protein = (row.get("Protein") or "").strip().upper()
            ligand = (row.get("Ligand") or "").strip()
            if not protein or not ligand:
                continue
            out_json = os.path.join(
                "results/boltz2_responses",
                f"{safe_name(protein)}_{safe_name(ligand)}.json",
            )
            if not os.path.exists(out_json):
                results.append(
                    {
                        "Protein": protein,
                        "Ligand": ligand,
                        "Affinity_Pred_Value": "",
                        "Affinity_pIC50": "",
                        "Affinity_Prob": "",
                        "Confidence": "",
                        "Endpoint": "",
                        "Status": "MISSING",
                    }
                )
                continue
            with open(out_json) as f:
                data = json.load(f)
            affinity_val, affinity_pic50, affinity_prob, confidence = parse_affinity(data)
            results.append(
                {
                    "Protein": protein,
                    "Ligand": ligand,
                    "Affinity_Pred_Value": affinity_val,
                    "Affinity_pIC50": affinity_pic50,
                    "Affinity_Prob": affinity_prob,
                    "Confidence": confidence,
                    "Endpoint": "",
                    "Status": "OK",
                }
            )
        out_csv = "results/boltz2_affinity_results.csv"
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
	                fieldnames=[
	                    "Protein",
	                    "Ligand",
	                    "Affinity_Pred_Value",
	                    "Affinity_pIC50",
	                    "Affinity_Prob",
	                    "Confidence",
	                    "Endpoint",
	                    "Status",
	                ],
	            )
            writer.writeheader()
            for r in results:
                writer.writerow(r)
        print(f"Wrote {out_csv}")
        return

    api_key = os.getenv("NVIDIA_API_KEY", "").strip()
    if not api_key:
        raise SystemExit(
            "Missing NVIDIA_API_KEY. Please set it in your environment, e.g.:\n"
            "  export NVIDIA_API_KEY='...'\n"
            "Do not hardcode API keys into scripts or packages."
        )

    results = []
    max_retries = int(os.getenv("BOLTZ2_MAX_RETRIES", "3"))
    backoff = float(os.getenv("BOLTZ2_RETRY_BACKOFF", "60"))
    for row in pairs:
        protein = (row.get("Protein") or "").strip().upper()
        ligand = (row.get("Ligand") or "").strip()
        if not protein or not ligand:
            continue
        if args.only_protein and protein != args.only_protein.upper():
            continue
        if args.only_ligand and ligand != args.only_ligand:
            continue

        uniprot_id = UNIPROT_IDS.get(protein)
        if not uniprot_id:
            print(f"[WARN] No UniProt ID for {protein}, skipping")
            continue

        smiles = smiles_map.get(ligand)
        if not smiles:
            smiles = normalized_smiles_map.get(normalize_ligand(ligand))
        if not smiles:
            print(f"[WARN] No SMILES for {ligand}, skipping")
            continue

        if protein == "MTOR":
            sequence = get_mtor_domain_sequence()
        else:
            sequence = fetch_uniprot_sequence(uniprot_id)

        response = None
        used_endpoint = None
        last_mode = ""
        last_endpoint = ""
        contact_text = (row.get("Contact_Residues") or "").strip()
        contacts = parse_contacts(contact_text) if args.use_pocket else []
        if not contacts:
            contacts = None

        out_json = os.path.join(
            "results/boltz2_responses",
            f"{safe_name(protein)}_{safe_name(ligand)}.json",
        )
        if args.resume and os.path.exists(out_json):
            results.append(
                {
                    "Protein": protein,
                    "Ligand": ligand,
                    "Affinity_Pred_Value": "",
                    "Affinity_Prob": "",
                    "Confidence": "",
                    "Endpoint": "",
                    "Status": "SKIP_EXISTS",
                }
            )
            continue

        auth_modes = [os.getenv("BOLTZ2_AUTH_MODE", "").strip().lower()]
        auth_modes = [m for m in auth_modes if m]
        if not auth_modes:
            auth_modes = ["bearer", "nvidia"]

        for endpoint in endpoints:
            try:
                for mode in auth_modes:
                    last_endpoint = endpoint
                    last_mode = mode
                    attempt = 0
                    response = None
                    while attempt < max_retries:
                        response = call_boltz2(
                            api_key,
                            endpoint,
                            sequence,
                            smiles,
                            contacts=contacts,
                            auth_mode=mode,
                        )
                        if response is not None and response.status_code == 429:
                            sleep_for = backoff * (attempt + 1)
                            print(f"[WARN] {protein}-{ligand} rate-limited. Sleeping {sleep_for:.0f}s.")
                            time.sleep(sleep_for)
                            attempt += 1
                            continue
                        break
                    if response is not None and response.status_code == 200:
                        used_endpoint = endpoint
                        break
                    if args.test and response is not None:
                        print(f"[TEST] endpoint={endpoint} mode={mode} status={response.status_code}")
                        print(response.text[:500])
                if response is not None and response.status_code == 200:
                    break
            except Exception as e:
                if os.getenv("BOLTZ2_VERBOSE"):
                    print(f"[DEBUG] exception={type(e).__name__} msg={e}")
                response = None
        if response is None or response.status_code != 200:
            status = response.status_code if response is not None else "ERR"
            if os.getenv("BOLTZ2_VERBOSE") and response is not None:
                print(f"[DEBUG] endpoint={last_endpoint} mode={last_mode} body={response.text[:500]}")
            print(f"[FAIL] {protein}-{ligand} status={status}")
            if response is not None:
                fail_path = os.path.join(
                    "results/boltz2_failures",
                    f"{safe_name(protein)}_{safe_name(ligand)}.txt",
                )
                with open(fail_path, "w") as f:
                    f.write(f"status={status}\n")
                    f.write(response.text)
            results.append(
                {
                    "Protein": protein,
                    "Ligand": ligand,
                    "Affinity_Pred_Value": "",
                    "Affinity_pIC50": "",
                    "Affinity_Prob": "",
                    "Confidence": "",
                    "Endpoint": used_endpoint or "",
                    "Status": status,
                }
            )
            continue

        data = response.json()
        affinity_val, affinity_pic50, affinity_prob, confidence = parse_affinity(data)

        with open(out_json, "w") as f:
            json.dump(data, f)

        results.append(
            {
                "Protein": protein,
                "Ligand": ligand,
                "Affinity_Pred_Value": affinity_val,
                "Affinity_pIC50": affinity_pic50,
                "Affinity_Prob": affinity_prob,
                "Confidence": confidence,
                "Endpoint": used_endpoint,
                "Status": response.status_code,
            }
        )
        if args.test:
            print("[TEST] Success. Exiting.")
            return
        time.sleep(args.sleep)

    out_csv = "results/boltz2_affinity_results.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Protein",
                "Ligand",
                "Affinity_Pred_Value",
                "Affinity_pIC50",
                "Affinity_Prob",
                "Confidence",
                "Endpoint",
                "Status",
            ],
        )
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
