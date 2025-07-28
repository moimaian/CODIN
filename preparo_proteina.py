import glob
import os
import subprocess
import shutil
from datetime import datetime
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QTextEdit

class PreparoProteina(QWidget):
    def __init__(self, codyn_dir: str, protein_dir: str, ligands_dir: str,
                 config_data: dict, status_prot: QTextEdit, run_folder: str, prot_path: str):
        super().__init__()
        self.codyn_dir = codyn_dir
        self.protein_dir = protein_dir
        self.ligands_dir = ligands_dir
        self.config_data = config_data
        self.status_prot = status_prot
        self.run_folder = run_folder
        self.prot_path = prot_path

    def log(self, message: str):
        """Adiciona mensagem no status_prot com timestamp e faz auto-scroll."""
        ts = datetime.now().strftime("%H:%M:%S")
        self.status_prot.append(f"[{ts}] {message}")
        self.status_prot.verticalScrollBar().setValue(
            self.status_prot.verticalScrollBar().maximum()
        )

    def edit_topol_top(self, run_folder, nome_ligante):
        topol_path = os.path.join(run_folder, "topol.top")
        if not os.path.exists(topol_path):
            self.log("Arquivo topol.top não encontrado para ajuste.")
            return

        with open(topol_path, "r") as f:
            lines = f.readlines()

        # 1. Inserir antes de [ moleculetype ]
        prm_block = [
            "; Include ligand parameters\n",
            f'#include "{nome_ligante}.prm"\n\n'
        ]

        # 2. Inserir antes de "; Include water topology"
        itp_block = [
            "; Include ligand topology\n",
            f'#include "{nome_ligante}.itp"\n\n'
        ]

        # Flags para controle de inserção
        new_lines = []
        prm_inserted = False
        itp_inserted = False

        for idx, line in enumerate(lines):
            # Inserir bloco prm antes de [ moleculetype ] (apenas na primeira ocorrência)
            if not prm_inserted and "[ moleculetype ]" in line:
                new_lines.extend(prm_block)
                prm_inserted = True
            # Inserir bloco itp antes de "; Include water topology" (apenas na primeira ocorrência)
            if not itp_inserted and "; Include water topology" in line:
                new_lines.extend(itp_block)
                itp_inserted = True
            new_lines.append(line)

        # Adicionar bloco ao final
        if not new_lines[-1].endswith("\n"):
            new_lines[-1] += "\n"
        new_lines.append(f"{nome_ligante}\t1\n")

        # Salva o arquivo ajustado
        with open(topol_path, "w") as f:
            f.writelines(new_lines)
        self.log(f"Arquivo topol.top ajustado para o ligante {nome_ligante}.")

    def start_preparation(self, nome_proteina: str, nome_ligante: str):
        try:
            self.log(f"\n####### COMPLEX PREPARATION START {nome_proteina} {nome_ligante}  #######")

            # Preparar diretórios
            self.log(f"Run directory ready: {self.run_folder}")

            # Copiar arquivos de TARGETS para a pasta de execução
            shutil.copy2(self.prot_path, self.run_folder)
            self.log("Copied protein into run folder.")

            # Carregar dados de configuração
            self.config_data = os.path.join(self.codyn_dir, ".form_data.txt")
            with open(self.config_data, 'r') as f:
                lines = f.readlines()
                ff = lines[0].strip()  # Campo de força
                solv = int(lines[6].strip())  # Modelo de Solvente
            self.log(f"Using Force Field: {ff}, solvent model: {solv}")
           
            # GROMACS preprocessing commands
            gmx = "/usr/local/gromacs/bin/gmx"
            self.log("Running gmx pdb2gmx")
            cmd = f"{gmx} pdb2gmx -f {nome_proteina}.pdb -o protein.gro"            
            subprocess.run(cmd, shell=True, cwd=self.run_folder, input=f"{ff}\n{solv}\n", text=True)

            # Unindo proteina e ligante
            self.log(f"Bringing {nome_ligante} and {nome_proteina} together in complex.gro")
            # Coletar o último número de resíduo da proteína
            resn_end = None
            with open(os.path.join(self.run_folder, "protein.gro")) as grof:
                for line in grof:
                    if len(line) >= 5 and line[:5].strip().isdigit():
                        resn = int(line[:5].strip())
                        resn_end = resn  # Vai sobrescrever até o último átomo
            self.log(f"Último número de resíduo da proteína: {resn_end}")
            resn_lig = resn_end + 1  # O ligante começa no próximo resíduo
            subprocess.run(f"{gmx} editconf -f {nome_ligante}.gro -o {nome_ligante}_new.gro -resnr {resn_lig}", shell=True, cwd=self.run_folder)

            # Remove header (2 linhas) e footer (última linha) do ligante
            self.log(f"Removing header and footer from {nome_ligante}_new.gro")
            subprocess.run(f"tail -n +3 {nome_ligante}_new.gro | head -n -1 > ligand_atoms.gro", shell=True, cwd=self.run_folder)

            # Remove footer do arquivo da proteína
            self.log(f"Removing header and footer from protein.gro")
            subprocess.run(f"tail -n +3 protein.gro | head -n -1  > protein_atoms.gro", shell=True, cwd=self.run_folder)

            # Junta os arquivos
            self.log("Combining protein and ligand atoms into complex_atoms.gro")
            subprocess.run(f"cat protein_atoms.gro ligand_atoms.gro > complex_atoms.gro", shell=True, cwd=self.run_folder)

            # Conta o número total de átomos
            self.log("Counting total number of atoms in complex_atoms.gro")
            result = subprocess.run("grep -v '^$' complex_atoms.gro | wc -l", shell=True, cwd=self.run_folder, capture_output=True, text=True)
            n_atoms = result.stdout.strip()
            self.log(f"Total number of atoms: {n_atoms}")

            # Caminhos dos arquivos
            complex_atoms_path = os.path.join(self.run_folder, "complex_atoms.gro")
            protein_gro_path = os.path.join(self.run_folder, "protein.gro")
            complex_gro_path = os.path.join(self.run_folder, "complex.gro")

            # Leia a última linha (box) do protein.gro
            with open(protein_gro_path, "r") as f:
                lines = f.readlines()
                box_line = lines[-1]

            # Leia os átomos do complex_atoms.gro
            with open(complex_atoms_path, "r") as f:
                atoms_data = f.read()

            # Escreva tudo em complex.gro
            with open(complex_gro_path, "w") as out:
                out.write(f"CODYN {nome_ligante} {nome_proteina}\n")
                out.write(f"{n_atoms}\n")
                out.write(atoms_data)
                if not atoms_data.endswith('\n'):
                    out.write('\n')
                out.write(box_line)

            self.edit_topol_top(self.run_folder, nome_ligante)

            # Finalização
            self.log("Complex preparation completed successfully.")

        except Exception as e:
            self.log(f"Error during protein prep: {e}")
