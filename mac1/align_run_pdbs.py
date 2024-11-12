from pymol import cmd
import glob

OUTDIRS = "OUTDIR_RUN_PDBS_ABF", "OUTDIR_RUN_PDBS_PNAS"

prefixes = []
for dirname in OUTDIRS:
    fnames = glob.glob(f"{dirname}/*.cif")

    prefixes += [cif.split(".")[0].upper() for cif in fnames]

# choose the reference
reference = f"{prefixes[0]}-monomer.pdb"

# align the others
for i_pref, pref in enumerate(prefixes[1:]):
    cmd.load(reference, "reference")
    print(f"aligning {pref} ({i_pref+1}/{len(prefixes[1:])})")
    obj_files = glob.glob(f"{pref}-*")
    for i_obj, obj_f in enumerate(obj_files):
        obj_name  = f"mobile{i_obj}"
        cmd.load(obj_f, obj_name)
        # use the segment name to distinguish objects
        cmd.alter(obj_name, f"segi={i_obj}")
    cmd.create("mobile", "all and not reference")
    cmd.align("mobile", "reference")
    
    for i_obj, obj_f in enumerate(obj_files):
        obj_name = f"mobile{i_obj}"
        cmd.create(obj_name, f"mobile and segi {i_obj} and not reference")
        cmd.save(obj_f, obj_name)
    cmd.delete("all")

print(f"Reference monomer: {reference}")

