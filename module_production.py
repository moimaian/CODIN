import glob
import os
import json
import subprocess
import shutil
from datetime import datetime
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import QWidget, QTextEdit

from PyQt5.QtCore import QThread, pyqtSignal
import subprocess

class MD_Thread(QThread):
    log_md_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)

    def __init__(self, gmx, run_folder, n_threads, t_ns_md):
        super().__init__()
        self.gmx = gmx
        self.run_folder = run_folder
        self.n_threads = n_threads
        self.t_ns_md = t_ns_md

    def run(self):
        # MD Production
        try:
            # Monta e roda o grompp para MD Production
            self.log_md_signal.emit("Running GROMPP for MD Production...")
            result = subprocess.run(
                f"{self.gmx} grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_'{self.t_ns_md}'.tpr",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_md_signal.emit(result.stdout.strip())
            if result.stderr:
                self.log_md_signal.emit("Erro: " + result.stderr.strip())

            # Roda MD Production (mdrun):
            self.log_md_signal.emit("Running MD Production...")
            result = subprocess.run(
                f"{self.gmx} mdrun -deffnm md_0_{self.t_ns_md} -ntmpi 1 -ntomp '{self.n_threads}' -nb gpu -bonded gpu -pme gpu -pmefft gpu -update gpu",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_md_signal.emit(result.stdout.strip())
            if result.stderr:
                self.log_md_signal.emit("Erro: " + result.stderr.strip())

            # Remove PBC (Periodic Boundary Conditions) from the trajectory
            self.log_md_signal.emit("Removing PBC from trajectory...")
            result = subprocess.run(
                f"{self.gmx} trjconv -s md_0_{self.t_ns_md}.tpr -f md_0_{self.t_ns_md}.xtc -o md_0_{self.t_ns_md}_noPBC.xtc -pbc mol -center",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_md_signal.emit(result.stdout.strip())            
            if result.stderr:
                self.log_md_signal.emit("Erro: " + result.stderr.strip())                
            self.finished_signal.emit(True)
            
        except Exception as e:
            self.log_md_signal.emit(f"Error during MD Production: {e}")
            self.finished_signal.emit(False)


class PrepMD(QWidget):
    log_md_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)
    def __init__(self, codyn_dir: str, config_data: dict, run_folder: str):
        super().__init__()
        self.codyn_dir = codyn_dir
        self.config_data = config_data
        self.run_folder = run_folder

    def log(self, message: str):
        """Adiciona mensagem no status_prot com timestamp e faz auto-scroll."""
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_md_signal.emit(f"[{ts}] {message}")

    def on_MD_finished(self, success):
        if success:
            self.log("Molecular Dynamics Finished Successfully!")

    def start_production(self, nome_proteina: str, nome_ligante: str):
        try:
            self.log(f"\n####### MD PRODUCTION START {nome_proteina} {nome_ligante} #######")

            # Carregar dados de configuração
            self.config_data = os.path.join(self.codyn_dir, ".form_data.txt")
            with open(self.config_data, 'r') as f:
                self.config = json.load(f)
                self.t_ns_md = self.config["production_time_ns"]
                self.n_threads = self.config["n_threads"]
            self.log(f"Using {self.n_threads} Threads in CPU for production")

            # GROMACS preprocessing commands
            gmx = "/usr/local/gromacs/bin/gmx"  

            # MD Production:
            self.log("Starting MD Production...")
            try:
                # Criando a thread de minimização
                self.MD_Thread = MD_Thread(gmx, self.run_folder, self.n_threads, self.t_ns_md)
                self.MD_Thread.log_md_signal.connect(self.log)
                self.MD_Thread.finished_signal.connect(self.on_MD_finished)
                self.MD_Thread.start()

            except Exception as e:
                self.log(f"Error during {nome_ligante} MD Production: {e}")

        except Exception as e:
            self.log(f"Error during MD Production: {e}")
