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

class NVT_Thread(QThread):
    log_eq_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)

    def __init__(self, gmx, run_folder, n_threads):
        super().__init__()
        self.gmx = gmx
        self.run_folder = run_folder
        self.n_threads = n_threads

    def run(self):
        # NVT equilibration
        try:
            # Monta e roda o grompp para NVT equilibration
            result = subprocess.run(
                f"{self.gmx} grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_eq_signal.emit(result.stdout.strip())
            if result.stderr:
                self.log_eq_signal.emit("Erro: " + result.stderr.strip())
            # Roda Equilibração NVT (mdrun):
            result = subprocess.run(
                f"{self.gmx} mdrun -deffnm nvt -ntmpi 1 -ntomp '{self.n_threads}' --nb gpu -bonded gpu -pme gpu -pmefft gpu -update gpu",
                shell=True, cwd=self.run_folder, capture_output=True, text=True
            )
            if result.stdout:
                self.log_eq_signal.emit(result.stdout.strip())
            
            if result.stderr:
                self.log_eq_signal.emit("Erro: " + result.stderr.strip())
            
            self.finished_signal.emit(True)
        
        except Exception as e:
            self.log_eq_signal.emit(f"Erro durante NVT equilibration: {e}")
            self.finished_signal.emit(False)


class PrepNVT(QWidget):
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

    def on_NVT_finished(self, success):
        if success:
            self.log("Minimização de energia finalizada com sucesso!")

            # Calculando a variaçao de temperatura:
            self.log("Calculating temperature variation...")
            gmx = "/usr/local/gromacs/bin/gmx"
            cmd = f"{gmx} energy -f nvt.edr -o temperature.xvg"
            result=subprocess.run(cmd, shell=True, cwd=self.run_folder, input="Temperature\n", text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

    def edit_topol_top2(self, run_folder, nome_ligante):
        topol_path = os.path.join(run_folder, "topol.top")
        if not os.path.exists(topol_path):
            self.log("Arquivo topol.top não encontrado para ajuste.")
            return

        with open(topol_path, "r") as f:
            lines = f.readlines()

        # Inserir position restraints
        pres_block = [
            "; Ligand position restraints\n"
            "#ifdef POSRES\n"
            '#include "posre_lig.itp"\n'
            "#endif\n\n"
        ]

        # Flags para controle de inserção
        new_lines = []
        pres_inserted = False

        for idx, line in enumerate(lines):
            # Inserir bloco prm antes de [ moleculetype ] (apenas na primeira ocorrência)
            if not pres_inserted and "; Include water topology" in line:
                new_lines.extend(pres_block)
                pres_inserted = True
            new_lines.append(line)

        # Salva o arquivo ajustado
        with open(topol_path, "w") as f:
            f.writelines(new_lines)
        self.log(f"Arquivo topol.top ajustado para o ligante {nome_ligante}.")

    def start_equilibration(self, nome_proteina: str, nome_ligante: str):
        try:
            self.log(f"\n####### EQUILIBRATION NVT START {nome_proteina} {nome_ligante} #######")

            # Carregar dados de configuração
            self.config_data = os.path.join(self.codyn_dir, "form_data.json")
            with open(self.config_data, 'r') as f:
                self.config = json.load(f)
                self.n_threads = self.config["n_threads"]
                p_ion = self.config["positive_ion"]
                n_ion = self.config["negative_ion"]
            self.log(f"Using {self.n_threads} Threads in CPU for equilibration")

            # GROMACS preprocessing commands
            gmx = "/usr/local/gromacs/bin/gmx"

            # Criando o arquivo de indices para o ligante:
            self.log("Creating a ligand Index file...")
            result=subprocess.run(f"echo '0 & ! a H*\nq' | {gmx} make_ndx -f {nome_ligante}.gro -o index_lig.ndx", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())
            
            # Criando o arquivo de restrições para o ligante:
            self.log("Creating a restrictions ligand file...")
            result=subprocess.run(f"echo '3' | {gmx} genrestr -f {nome_ligante}.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())

            self.edit_topol_top2(self.run_folder, nome_ligante)
            
            # Criando os acoplamentos de temperatura:
            self.log("Creating tc-groups...")
            result=subprocess.run(f"echo '\"Protein\" | \"{nome_ligante}\"\n\"{p_ion}\" | \"{n_ion}\" | \"SOL\"\nq' | {gmx} make_ndx -f em.gro -o index.ndx", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            if result.stdout:
                self.log(result.stdout.strip())
            if result.stderr:
                self.log("Erro: " + result.stderr.strip())            
            

            # NVT Equilibration:
            self.log("Starting NVT...")
            try:
                # Criando a thread de minimização
                self.nvt_thread = NVT_Thread(gmx, self.run_folder, self.n_threads)
                self.nvt_thread.log_eq_signal.connect(self.log)
                self.nvt_thread.finished_signal.connect(self.on_NVT_finished)
                self.nvt_thread.start()

            except Exception as e:
                self.log(f"Error during {nome_ligante} NVT Equilibration: {e}")            

        except Exception as e:
            self.log(f"Error during Equilibration: {e}")
