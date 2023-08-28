# Read in a Dr. T output file and just pull out the most prevelent structure present
# at the end of the simulated transcription reaction

from pathlib import Path
import re
import pandas as pd
import numpy as np

TRIGGER_LINE = '# Distribution of structures at the end:'


def main():

    logfile = []
    for each_file in Path(snakemake.input[0]).iterdir():
        if each_file.suffix == '.log':
            logfile = each_file
            break

    with open(str(logfile)) as handle:
        struct_stats = [[], []]  # energy and occupancy
        read_entries = False
        line = handle.readline()
        while line:
            if read_entries == False:
                if line.strip() == TRIGGER_LINE:
                    # skip one line
                    handle.readline()
                    read_entries = True
            else:
                matches = re.findall(r'\d+\.\d+', line)
                struct_stats[0].append(float(matches[0]))
                struct_stats[1].append(float(matches[1]))
            line = handle.readline()
    
    struct_stats = np.array(struct_stats)
    df_stats = pd.DataFrame({
        'energy': struct_stats[0],
        'occupany': struct_stats[1]
    })
    df_stats = df_stats.loc[df_stats.occupany > 0.01]
    df_stats = df_stats.sort_values(by=['occupany'])
    
    structre_energy = df_stats.iloc[-1].energy

    with open(snakemake.output[0], 'w') as write_path:
        write_path.write(str(structre_energy))


if __name__ == '__main__':
    main()
            




