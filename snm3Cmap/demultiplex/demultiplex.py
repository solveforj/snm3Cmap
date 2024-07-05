from pathlib import Path
import os
import yaml

class PrepareDemultiplex:

    def snm3Cseq(self):

        barcodes = self.config_dict["setup"]["barcodes"]
        plate_info = self.config_dict["setup"]["plate_info"]
        output_directory = self.config_dict["setup"]["output_directory"]

        Path(output_directory).mkdir(parents=True, exist_ok=True)
        
        if self.mode == "snm3Cseq":
            smk_path = 'demultiplex_snm3Cseq.Snakefile'
        elif self.mode == "scalemethyl":
            smk_path = 'demultiplex_scalemethyl.Snakefile'
            
        smk_path = os.path.join(Path(__file__).parent.resolve(), "smk", smk_path)
    
        with open(smk_path) as f:
            snake_template = f.read()
                    
        with open(plate_info) as f:
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == "#":
                    continue
                line = line.split()
                plate, fastq_directory = line[0], line[1]
                plate_run_directory = os.path.join(output_directory, plate)
                Path(plate_run_directory).mkdir(parents=True, exist_ok=True)
    
                plate_config = {}
                for mate in ["R1", "R2"]: 
                    r_fn = f"{fastq_directory}*{plate}[-_]*{mate}*fastq.gz"
                    plate_config[mate] = r_fn
    
                plate_config["out_dir"] = plate_run_directory
                plate_config["barcodes"] = barcodes
                plate_config["plate"] = plate
    
                yaml_directory = f"{plate_run_directory}/config.yaml"
                
                with open(yaml_directory, "w") as outfile:
                    yaml.dump(plate_config, outfile)
    
                with open(f"{plate_run_directory}/demultiplex.smk", 'w') as outfile:
                    outfile.write(snake_template)
                    
                with open(f"{plate_run_directory}/demultiplex_cmd.txt", 'w') as outfile:
                    cmd = f"snakemake -d {plate_run_directory} "
                    cmd += f"--snakefile {plate_run_directory}/demultiplex.smk "
                    cmd += f" -c {self.jobs} "
                    cmd += f" {self.nolock} {self.rerun_incomplete} "
                    outfile.write(cmd + '\n')

    def __init__(self,
                 config,
                 jobs=2,
                 nolock=False,
                 rerun_incomplete=False
                ):
        
        with open(config) as f:
            self.config_dict = yaml.safe_load(f)

        self.mode = self.config_dict["setup"]["mode"]
        
        self.jobs = jobs
        
        self.nolock = "--nolock" if nolock else ""
        self.rerun_incomplete = "--rerun-incomplete" if rerun_incomplete else ""

        if "snm3Cseq" == self.mode:
            self.snm3Cseq()
        elif "scalemethyl" == self.mode:
            self.snm3Cseq()
        else:
            print("Mode not recognized!")