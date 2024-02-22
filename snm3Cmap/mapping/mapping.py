from pathlib import Path
from glob import glob
import os
import yaml

def prepare_mapping(plate_info,
                    demultiplex_directory, 
                    barcodes,
                    reference_genome,
                    jobs=2,
                    nolock=False,
                    rerun_incomplete=False
                   ):

    smk_path = 'mapping.Snakefile'
    
    with Path(__file__).with_name(smk_path).open('r') as f:
        snake_template = f.read()
    
    cell_ids = []
    with open(barcodes) as f:
        for line in f:
            line = line.strip()
            if line[0] == ">":
                cell_ids.append(line[1:])

    biscuit_threads = jobs // 2 if jobs > 1 else 1
    
    params = {
        "biscuit_threads" : f'{biscuit_threads}',
        "reference_path" : f'"{reference_genome}"'
    }

    nolock = "--nolock" if nolock else ""
    rerun_incomplete = "--rerun-incomplete" if rerun_incomplete else ""
    
    with open(plate_info) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            line = line.split("\t")
            fastq_directory, plate = line[0], line[1]
            
            plate_run_directory = os.path.join(demultiplex_directory, plate)

            with open(f"{plate_run_directory}/mapping_scripts.txt", 'w') as scripts:
                
                for cell in cell_ids:
                    cell_local_directory = f"{plate}_{cell}"
                    cell_run_directory = os.path.join(plate_run_directory, 
                                                      cell_local_directory)
                    
                    snakemake_cmd = f"{cell_run_directory}/mapping_cmd.txt"
                    scripts.write(snakemake_cmd + "\n")
                    
                    with open(snakemake_cmd, 'w') as f:
                        cmd = f"snakemake -d {cell_run_directory} "
                        cmd += f"--snakefile {plate_run_directory}/mapping.smk -c {jobs} "
                        cmd += f"--config cell={cell} plate={plate} {nolock} {rerun_incomplete} "
                        f.write(cmd + '\n')

            params_write = "\n".join([f"{k} = {v}" for k, v in params.items()])
            params_write += "\n"

            with open(f"{plate_run_directory}/mapping.smk", 'w') as f:
                f.write(params_write + snake_template)
                
