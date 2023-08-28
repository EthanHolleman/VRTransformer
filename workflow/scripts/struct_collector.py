# Given output directory

from pathlib import Path
import pandas as pd


def collect_structure_energies(dir: Path):
    energy_files = []

    def walk_dir(dir):
        if dir.is_dir():
            for each_path in dir.iterdir():
                if each_path.is_dir():
                    walk_dir(each_path)
                elif each_path.suffix == ".energy":
                    energy_files.append(each_path)
        elif each_path.suffix == ".energy":
            energy_files.append(each_path)

    walk_dir(dir)

    return energy_files


def make_energy_dicts(energy_files: list):
    energy_dicts = []

    def assign_energy_to_permutation(energy_file: Path):
        # parent directory of the energy file should always be the permutation number
        # based on snakemake rules
        perm = Path(energy_file).parent.parts[-1]

        with open(str(energy_file)) as handle:
            energy = float(handle.readline().strip())

        return {
            "permutation": perm,
            'perm_id': perm.split('-')[-1],
            "energy": energy}

    for each_file in energy_files:
        energy_dicts.append(assign_energy_to_permutation(each_file))

    return energy_dicts


def add_fasta_path_and_plasmid_name_to_energy_dicts(
    plasmid_name: str, energy_dicts: list, permutation_fa_list: list
):
    fa_lookup = {Path(p).stem: p for p in permutation_fa_list}

    def determine_perm_fasta(perm_dict: dict):
        fasta_path = fa_lookup[perm_dict["permutation"]]
        perm_dict["fasta_path"] = str(fasta_path)

        return perm_dict

    for i in range(len(energy_dicts)):
        energy_dicts[i] = determine_perm_fasta(energy_dicts[i])
        energy_dicts[i]['plasmid'] = plasmid_name

    return energy_dicts


def main():
    energy_dicts = make_energy_dicts(snakemake.input['energy_files'])
    energy_dicts = add_fasta_path_and_plasmid_name_to_energy_dicts(
        snakemake.params['plasmid'], energy_dicts, snakemake.input["plasmids"]
    )

    pd.DataFrame(energy_dicts).to_csv(snakemake.output[0], sep="\t", index=False)


if __name__ == "__main__":
    main()
