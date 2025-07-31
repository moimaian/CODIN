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

class NPT_Thread(QThread):
    log_eq_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)

    def __init__(self, gmx, run_folder, n_threads):
        super().__init__()
        self.gmx = gmx
        self.run_folder = run_folder
        self.n_threads = n_threads

    def run(self):
        # NPT equilibration
        try:
            # Monta e roda o grompp para NPT equilibration
            result = subprocess.run(
                f"{self.gmx} grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_eq_signal.emit(result.stdout.strip())
            if result.stderr:
                self.log_eq_signal.emit("Erro: " + result.stderr.strip())

            # Roda Equilibração NPT (mdrun):
            result = subprocess.run(
                f"{self.gmx} mdrun -deffnm npt -ntmpi 1 -ntomp '{self.n_threads}' -nb gpu -bonded gpu -pme gpu -pmefft gpu -update gpu",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_eq_signal.emit(result.stdout.strip())
            
            if result.stderr:
                self.log_eq_signal.emit("Erro: " + result.stderr.strip())            
            self.finished_signal.emit(True)
        except Exception as e:
            self.log_eq_signal.emit(f"Erro durante NPT equilibration: {e}")
            self.finished_signal.emit(False)


class PrepNPT(QWidget):
    log_eq_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)
    def __init__(self, codyn_dir: str, config_data: dict, run_folder: str):
        super().__init__()
        self.codyn_dir = codyn_dir
        self.config_data = config_data
        self.run_folder = run_folder

    def log(self, message: str):
        """Adiciona mensagem no status_prot com timestamp e faz auto-scroll."""
        ts = datetime.now().strftime("%H:%M:%S")
        self.log_eq_signal.emit(f"[{ts}] {message}")

    def on_NPT_finished(self, success):
        if success:
            self.log("Minimização de energia finalizada com sucesso!")

            # Calculando a variaçao de pressão:
            self.log("Calculating pressure variation...")
            gmx = "/usr/local/gromacs/bin/gmx"
            cmd = f"{gmx} energy -f npt.edr -o pressure.xvg"
            result=subprocess.run(cmd, shell=True, cwd=self.run_folder, input="Pressure\n", text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

            # Exibindo o gráfico da variação de pressão:
            self.log("Displaying pressure variation graph...")
            if shutil.which("xmgrace"):
                subprocess.run("xmgrace pressure.xvg -autoscale xy", shell=True, cwd=self.run_folder)

            # Calculando a variaçao de densidade:
            self.log("Calculating density variation...")
            gmx = "/usr/local/gromacs/bin/gmx"
            cmd = f"{gmx} energy -f npt.edr -o density.xvg"
            result=subprocess.run(cmd, shell=True, cwd=self.run_folder, input="Density\n", text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

            # Exibindo o gráfico da variação de densidade:
            self.log("Displaying density variation graph...")
            if shutil.which("xmgrace"):
                subprocess.run("xmgrace density.xvg -autoscale xy", shell=True, cwd=self.run_folder)
            self.log("\nNPT Equilibration completed successfully!\n")
        else:
            self.log("Error during execution of NPT Equilibration.")


    def start_equilibration(self, nome_proteina: str, nome_ligante: str):
        try:
            self.log(f"\n####### EQUILIBRATION NPT START {nome_proteina} {nome_ligante} #######")

            # Carregar dados de configuração
            self.config_data = os.path.join(self.codyn_dir, "form_data.json")
            with open(self.config_data, 'r') as f:
                self.config = json.load(f)
                self.n_threads = self.config["n_threads"]
            self.log(f"Using {self.n_threads} Threads in CPU for equilibration")

            # GROMACS preprocessing commands
            gmx = "/usr/local/gromacs/bin/gmx"  
            
            # NPT Equilibration:
            self.log("Starting NPT...")
            try:
                # Criando a thread de minimização
                self.npt_thread = NPT_Thread(gmx, self.run_folder, self.n_threads)
                self.npt_thread.log_eq_signal.connect(self.log)
                self.npt_thread.finished_signal.connect(self.on_NPT_finished)
                self.npt_thread.start()

            except Exception as e:
                self.log(f"Error during {nome_ligante} NPT Equilibration: {e}")

        except Exception as e:
            self.log(f"Error during Equilibration: {e}")
