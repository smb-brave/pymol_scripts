
from pymol import cmd

try:
    import numpy as np
    import pandas
    from pypdb import describe_chemical
except ImportError:
    import sys
    install_cmd = f"{sys.executable} -m pip install numpy pandas pypdb"
    raise ImportError(f"Install missing dependencies using `{install_cmd}`")


def _B_and_q(obj):
    """get average B factor and occupancy for pymol object"""
    Bs = []
    cmd.iterate(obj, lambda atom: Bs.append(atom.b))
    ave_B = np.mean(Bs)

    qs = []
    cmd.iterate(obj, lambda atom: qs.append(atom.q))
    qs = set(qs)
    assert len(qs)==1
    q = list(qs)[0]
    return ave_B, q


def run_pdbs(df, outdir, csv_f):
    """
    df: pandas.DataFrame with a `pdb_code` column
    outdir, str outputfolder
    csv_f, str output filename
    """

    os.makedirs(outdir, exist_ok=True)

    all_monomer_files = []
    all_ligand_records = []
    check_me = []
    for i_p, p in enumerate(df.pdb_code):
        monomer_files = []

        print(f"Beginning parse of {p} ({i_p+1}/{len(df)})")
        # download the CIF file for pdb code `p`
        cmd.fetch(p, path=outdir)
        
        # create objects representing monomers and ligands separately
        cmd.create("monomers", "polymer.protein")
        cmd.create("ligs", "organic and not solvent")

        # most mac1 pdbs are of the dimer, occasionally there is only one monomer in the pdb
        # so expect chain labels A,B
        chain_names = []
        cmd.iterate("monomers", lambda atom: chain_names.append(atom.chain))
        chain_names = set(chain_names)
        assert "A" in chain_names or "B" in chain_names and len(chain_names) in [1,2]
        is_dimer =  len(chain_names) == 2 
        if not is_dimer:
            assert chain_names == {"A"}

        # ligands have different residue names
        lig_names = []
        cmd.iterate("ligs", lambda atom: lig_names.append(atom.resn))
        lig_names = set(lig_names)

        # list of SDF files, Bfacs, and occupancy for the ligands (multiple conformers)
        this_pdb_lig_records = []

        # loop over ligand residue names `resn`
        for resn in lig_names:
            # skip the dms molecule  (dimethyl sulfoxide)
            if resn=="DMS":
                continue
            # create object for this ligand `resn`
            cmd.create(f"lig{resn}", f"ligs and resn {resn}")

            # list occupancies for this ligand `resn`, occasionally there are > 1
            occ = []
            cmd.iterate(f"lig{resn}", lambda atom: occ.append(atom.q))
            occ = [round(q, 2) for q in set(occ)]
            print(p, resn, occ)
            
            # segments representing this ligand `resn`
            lig_segments = []
            cmd.iterate(f"lig{resn}", lambda atom: lig_segments.append(atom.segi))
            lig_segments = set(lig_segments)

            # sometimes multiple conformors have same segment ID...  in that case , use occupancy to iterate over conformers... 
            # note , if occupancy is somehow 0.5, then this might skip one conformer!
            if len(occ) > len(lig_segments):
                lig_iter = enumerate(occ)
                use_q = True
            else:
                lig_iter = enumerate(lig_segments)
                use_q = False

            for i_lig, val in lig_iter:
                if use_q:
                    lig_sel = f"organic and not solvent and q={val} and resn {resn}"
                else:
                    lig_sel = f"organic and not solvent and segi {val} and resn {resn}"

                lig_obj = f"lig{resn}-{i_lig}" 
                cmd.create(lig_obj, lig_sel) 
                lig_f = f"{outdir}/{p}-lig{resn}-{i_lig}.sdf"
                cmd.save(lig_f, lig_obj)
                # get b factor and q
                B, q = _B_and_q(lig_obj)
                lig_record = f"{lig_f}:{B:.2f}:{q:.2f}"
                this_pdb_lig_records.append(lig_record)
                if use_q:
                    assert np.allclose(val, q)
                    if q == 0.5:
                        check_me.append(p)

            
            # extract the monomer closest to the deposited ligands for docking and posebusters validation
            if is_dimer:
                # if a dimer, find the monomer nearest to the ligands
                com_A = np.array(cmd.centerofmass("chain A and polymer.protein"))
                com_B = np.array(cmd.centerofmass("chain B and polymer.protein"))
                com_L = np.array(cmd.centerofmass("organic and not solvent"))
                dist_A = np.sqrt(sum((com_A-com_L)**2))
                dist_B = np.sqrt(sum((com_B-com_L)**2))
                dists = {"A":dist_A, "B": dist_B}
                closest_chain = "A" if dist_A < dist_B else "B"
            else:
                closest_chain = "A"
            cmd.create("monomer", f"polymer.protein and not solvent and not organic and chain {closest_chain}")
            monomer_f = f"{outdir}/{p}-monomer.pdb"
            cmd.save(monomer_f, "monomer")

        all_monomer_files.append(monomer_f)
        all_ligand_records.append (";".join(this_pdb_lig_records))

        cmd.delete('all')

    df["monomer_pdbs"] = all_monomer_files 
    df["ligand_records"] = all_ligand_records
    df = df.assign(ligand_records=df['ligand_records'].str.split(';')).explode('ligand_records').reset_index(drop=True)
    df[["ligand_sdfs","ligand_Bfac","ligand_occ"]] = df.ligand_records.str.split(":", n=2, expand=True)
    lig_code = [l.split("/")[1].split("-lig")[1].split("-")[0] for l in df.ligand_sdfs.values]
    df['lig_code'] = lig_code
    df.to_csv(csv_f, sep="\t", index=False)
    return check_me


def _add_rcsb_smiles(csv_f):
    """
    csv_f: tab separated database with a column `lig_code`, 3-digit code representing
    ligands in the PDB

    add smiles and sterosmiles strings from the PDB (sometimes they are different due
    to sterochemistry), then resaves the `csv_f` file 
    this uses the describe_chemical method from pypdb
    """
    df = pandas.read_csv(csv_f, sep="\t")

    smis, stereo_smis = [],[]
    for i_lc, lc in enumerate(df.lig_code):
        lig_dat = describe_chemical(lc)
        descr = lig_dat['rcsb_chem_comp_descriptor']
        smi = descr['smiles']
        stereo_smi = descr['smilesstereo']
        smis.append(smi)
        stereo_smis.append(stereo_smi)
        print(i_lc+1, len(df), stereo_smi)

    df["rcsb_smiles"] = smis
    df["rcsb_stereosmiles"] = stereo_smis
    df.to_csv(csv_f, sep="\t", index=False)

# read the schuller et al database
df_abs = pandas.read_excel("abf8711_data_file_s1.xlsx", sheet_name="11. fragment_hits_classified", skiprows=1)
df_abs["pdb_code"] = df_abs.PDB
outdir = "OUTDIR_RUN_PDBS_ABF"
csv_f = "RUN_PDBS_ABF.csv"
check_me = run_pdbs(df_abs, outdir, csv_f)
_add_rcsb_smiles(csv_f)

# read the pnas database 
df = pandas.read_excel("pnas.2212931120.sd01.xlsx", sheet_name="Compound Summary")
df = df.loc[df.PDB.notnull()]
gb = df.groupby("PDB")
c = gb.count()
df = [gb.get_group(p) for p in c.index[ c.Smiles == 1] ]
df = pandas.concat(df)
pdbs = [ p.strip().split()[0] for p in df.PDB]
df["pdb_code"] = pdbs
outdir = "OUTDIR_RUN_PDBS_PNAS"
csv_f = "RUN_PDBS_PNAS.csv"
check_me += run_pdbs(df, outdir, csv_f)
_add_rcsb_smiles(csv_f)

print("Double check these PDBs for multple conformers:")

