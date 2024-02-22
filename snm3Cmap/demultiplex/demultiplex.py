from pathlib import Path
import os
import yaml

def prepare_demultiplex(plate_info,
                        output_directory,
                        barcodes,
                        jobs=2,
                        nolock=False,
                        rerun_incomplete=False
                       ):

    Path(output_directory).mkdir(parents=True, exist_ok=True)
    
    with Path(__file__).with_name('demultiplex.Snakefile').open('r') as f:
        snake_template = f.read()
        
    nolock = "--nolock" if nolock else ""
    rerun_incomplete = "--rerun-incomplete" if rerun_incomplete else ""

    with open(plate_info) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            line = line.split("\t")
            fastq_directory, plate = line[0], line[1]
            plate_run_directory = os.path.join(output_directory, plate)
            Path(plate_run_directory).mkdir(parents=True, exist_ok=True)

            plate_config = {}
            for mate in ["R1", "R2"]: 
                r_fn = f"{fastq_directory}*{plate}*{mate}*fastq.gz"
                plate_config[mate] = r_fn

            plate_config["out_dir"] = plate_run_directory
            plate_config["barcodes"] = barcodes
            plate_config["plate"] = plate

            yaml_directory = f"{plate_run_directory}/config.yaml"
            with open(yaml_directory, "w") as outfile:
                yaml.dump(plate_config, outfile)

            with open(f"{plate_run_directory}/demultiplex.smk", 'w') as f:
                f.write(snake_template)
                
            with open(f"{plate_run_directory}/demultiplex_cmd.txt", 'w') as f:
                cmd = f"snakemake -d {plate_run_directory} "
                cmd += f"--snakefile {plate_run_directory}/demultiplex.smk "
                cmd += f" -c {jobs} "
                cmd += f" {nolock} {rerun_incomplete} "
                f.write(cmd + '\n')
