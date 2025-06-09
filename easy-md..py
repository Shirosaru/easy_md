import os
import requests
import subprocess
from easy_md.main.quickrun import quickrun
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def download_pdb(pdb_id, out_file="protein.pdb"):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"📥 Downloading PDB from {url}")
    r = requests.get(url)
    if r.status_code != 200:
        raise Exception(f"❌ Failed to download PDB: {pdb_id}")
    with open(out_file, "w") as f:
        f.write(r.text)
    print(f"✅ Saved as {out_file}")
    return out_file

def fix_protein(in_file, out_file="protein_fixed.pdb"):
    print(f"🔧 Fixing PDB file {in_file} with pdbfixer...")
    fixer = PDBFixer(filename=in_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(out_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"✅ Fixed PDB saved as {out_file}")
    return out_file

def remove_hetatm_from_pdb(pdb_file, resnames_to_remove, out_file="protein_cleaned.pdb"):
    resnames_to_remove = {r.upper() for r in resnames_to_remove}
    with open(pdb_file, "r") as fin, open(out_file, "w") as fout:
        for line in fin:
            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                if resname in resnames_to_remove:
                    continue
            fout.write(line)
    print(f"🧹 Removed residues: {', '.join(resnames_to_remove)}")
    return out_file

def run_md_simulation(protein_file, ligand_file=None, nsteps=1000):
    print(f"🚀 Running MD simulation...")
    try:
        quickrun(
            protein_file=protein_file,
            ligand_file=ligand_file,
            nsteps=nsteps
        )
        print("✅ MD completed.")
    except Exception as e:
        print(f"❌ MD failed: {e}")

if __name__ == "__main__":
    pdb_id = "6lkd"
    raw_pdb = download_pdb(pdb_id)
    fixed_pdb = fix_protein(raw_pdb)

    # Remove unwanted hetero groups (FAD, EGO, SO4)
    cleaned_pdb = remove_hetatm_from_pdb(fixed_pdb, {"FAD", "EGO", "SO4"})

    # Run simulation without ligands or cofactors
    run_md_simulation(cleaned_pdb, ligand_file=None, nsteps=1000)
