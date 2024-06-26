from pathlib import Path
from glob import glob
import os
import yaml

class PrepareMapping:

    def snm3Cseq(self):

        barcodes = self.config_dict["setup"]["barcodes"]
        plate_info = self.config_dict["setup"]["plate_info"]
        output_directory = self.config_dict["setup"]["output_directory"]

        Path(output_directory).mkdir(parents=True, exist_ok=True)
        
        if self.mode == "snm3Cseq":
            smk_path = 'mapping_snm3Cseq.Snakefile'
            
        smk_path = os.path.join(Path(__file__).parent.resolve(), "smk", smk_path)
    
        with open(smk_path) as f:
            snake_template = f.read()
                    
        cell_ids = []
        with open(barcodes) as f:
            for line in f:
                line = line.strip()
                if line[0] == ">":
                    cell_ids.append(line[1:])
        
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
    
                with open(f"{plate_run_directory}/mapping_scripts.txt", 'w') as scripts:
                    
                    for cell in cell_ids:
                        cell_local_directory = f"{plate}_{cell}"
                        cell_run_directory = os.path.join(plate_run_directory, 
                                                          cell_local_directory)
                        
                        snakemake_cmd = f"{cell_run_directory}/mapping_cmd.txt"
                        scripts.write(snakemake_cmd + "\n")
                        
                        with open(snakemake_cmd, 'w') as f:
                            cmd = f"snakemake -d {cell_run_directory} "
                            cmd += f"--snakefile {plate_run_directory}/mapping.smk -c {self.jobs} "
                            cmd += f"--config cell={cell} plate={plate} {self.nolock} {self.rerun_incomplete} "
                            f.write(cmd + '\n')
            
                params_write = "\n".join([f"{k} = {v}" for k, v in self.params.items()])
                params_write += "\n"
            
                with open(f"{plate_run_directory}/mapping.smk", 'w') as f:
                    f.write(params_write + snake_template)
    
    def demultiplexed(self):

        if self.mode == "contacts":
            smk_path = 'mapping_contacts.Snakefile'
        elif self.mode == "contacts_bisulfite":
            smk_path = 'mapping_contacts_bisulfite.Snakefile'
        else:
            raise ValueError('Valid mode was not provided!')

        smk_path = os.path.join(Path(__file__).parent.resolve(), "smk", smk_path)
        
        output_directory = self.config_dict["setup"]["output_directory"]
        fastq_info = self.config_dict["setup"]["fastq_info"]

        Path(output_directory).mkdir(parents=True, exist_ok=True)
        
        with open(smk_path) as f:
            snake_template = f.read()
            
        with open(f"{output_directory}/mapping_scripts.txt", 'w') as scripts, \
            open(fastq_info) as f:
    
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == "#":
                    continue
                line = line.split()
                exp, r1_fastq, r2_fastq = line[0], line[1], line[2]
                
                exp_run_directory = os.path.join(output_directory, exp)
                exp_run_directory_path = Path(exp_run_directory)
                exp_run_directory_path.mkdir(parents=True, exist_ok=True)
                        
                snakemake_cmd = f"{exp_run_directory}/mapping_cmd.txt"
                scripts.write(snakemake_cmd + "\n")
                        
                with open(snakemake_cmd, 'w') as f:
                    cmd = f"snakemake -d {exp_run_directory} "
                    cmd += f"--snakefile {output_directory}/mapping.smk -c {self.jobs} "
                    cmd += f"--config r1={r1_fastq} r2={r2_fastq} exp={exp} {self.nolock} {self.rerun_incomplete} "
                    f.write(cmd + '\n')

        params_write = "\n".join([f"{k} = {v}" for k, v in self.params.items()])
        params_write += "\n"
            
        with open(f"{output_directory}/mapping.smk", 'w') as f:
            f.write(params_write + snake_template)

    def __init__(self,
                 config,
                 jobs=2,
                 nolock=False,
                 rerun_incomplete=False
                ):
        
        with open(config) as f:
            self.config_dict = yaml.safe_load(f)

        self.mode = self.config_dict["setup"]["mode"]
        
        mapping_threads = jobs // 2 if jobs > 1 else 1

        self.jobs = jobs
        
        self.params = {
            "mapping_threads" : f'{mapping_threads}',
            "parameters_path" : f'"{config}"',
        }
    
        self.nolock = "--nolock" if nolock else ""
        self.rerun_incomplete = "--rerun-incomplete" if rerun_incomplete else ""

        if "snm3Cseq" in self.mode:
            
            self.snm3Cseq()

        else:

            self.demultiplexed()

                
