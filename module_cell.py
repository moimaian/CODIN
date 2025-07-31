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

class MinimizationThread(QThread):
    log_cell_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)

    def __init__(self, gmx, run_folder):
        super().__init__()
        self.gmx = gmx
        self.run_folder = run_folder

    def run(self):
        # Minimização de energia
        try:
            # Monta e roda o grompp para energia minimization
            result = subprocess.run(
                f"{self.gmx} grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_cell_signal.emit(result.stdout.strip())
            if result.stderr:
                self.log_cell_signal.emit("Erro: " + result.stderr.strip())
                
            # Roda minimização (mdrun)
            result = subprocess.run(
                f"{self.gmx} mdrun -v -deffnm em",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_cell_signal.emit(result.stdout.strip())
            
            if result.stderr:
                self.log_cell_signal.emit("Erro: " + result.stderr.strip())
                
            self.finished_signal.emit(True)
            
        except Exception as e:
            self.log_cell_signal.emit(f"Erro durante minimização: {e}")
            self.finished_signal.emit(False)


class PreparoCell(QWidget):
    log_cell_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)
    def __init__(self, codyn_dir: str, config_data: dict, run_folder: str):
        super().__init__()
        self.codyn_dir = codyn_dir
        self.config_data = config_data
        self.run_folder = run_folder

    def log(self, message: str):
        """Adiciona mensagem no status_prot com timestamp e faz auto-scroll."""
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_cell_signal.emit(f"[{ts}] {message}")

    def on_minimization_finished(self, success):
        if success:
            self.log("\nEnergy Minimization completed successfully!\n")

            # Calculando a energia potencial:
            self.log("Calculating potential energy...")
            gmx = "/usr/local/gromacs/bin/gmx"
            cmd = f"{gmx} energy -f em.edr -o potential.xvg"
            result=subprocess.run(cmd, shell=True, cwd=self.run_folder, input="Potential\n", text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

            # Exibindo o gráfico da energia potencial:
            self.log("Displaying potential energy graph...")
            if shutil.which("xmgrace"):
                subprocess.run("xmgrace potential.xvg -autoscale xy", shell=True, cwd=self.run_folder)            

        else:
            self.log("Error during execution of energy minimization.")  

    def start_preparation(self, nome_proteina: str, nome_ligante: str):
        try:
            self.log(f"\n####### CELL PREPARATION START {nome_proteina} {nome_ligante} #######")

            # Carregar dados de configuração
            self.config_data = os.path.join(self.codyn_dir, "form_data.json")
            with open(self.config_data, 'r') as f:
                self.config = json.load(f)
                conc_ions = self.config["ion_concentration_M"]
                p_ion = self.config["positive_ion"]
                n_ion = self.config["negative_ion"]
            self.log(f"Using ion conc: {conc_ions} M, positive ion: {p_ion}, negative ion: {n_ion}")
           
            # GROMACS preprocessing commands
            gmx = "/usr/local/gromacs/bin/gmx"

            # Criando a célula unitária:
            self.log("Creating unit cell...")
            result=subprocess.run(f"{gmx} editconf -f complex.gro -o newbox.gro -bt cubic -c -d 1.0", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())
            
            # Preenchendo a célula unitária com solvente:
            self.log("Filling unit cell with solvent...")
            result=subprocess.run(f"{gmx} solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

            # Preparando a topologia dos íons:
            self.log("Preparing ion topology...")
            result=subprocess.run(f"{gmx} grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

            # Adding ions:
            self.log(f"Adding ions {p_ion} and {n_ion} ...")
            subprocess.run(f"echo 'SOL' | {gmx} genion -s ions.tpr -o solv_ions.gro -p topol.top -pname {p_ion} -nname {n_ion} -conc {conc_ions} -neutral", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())
            self.log("\nCell preparation completed successfully!\n")

            # Energy Minimization:
            self.log("Starting energy minimization...")
            try:
                # Criando a thread de minimização
                self.min_thread = MinimizationThread(gmx, self.run_folder)
                self.min_thread.log_cell_signal.connect(self.log)
                self.min_thread.finished_signal.connect(self.on_minimization_finished)
                self.min_thread.start()

            except Exception as e:
                self.log(f"Error during starting energy minimization: {e}")            

        except Exception as e:
            self.log(f"Error during cell preparation: {e}")
