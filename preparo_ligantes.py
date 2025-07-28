import os
import subprocess
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QLabel, QMessageBox, QSizePolicy
)
from PyQt5.QtGui import QPixmap
from datetime import datetime
import shutil
import numpy as np
from PyQt5.QtWidgets import QTextEdit
import subprocess
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    import networkx as nx
except ImportError:
    print("networkx não encontrado. Instalando...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "networkx"])
    import networkx as nx

class PreparoLigantes(QWidget):
    def __init__(self, codyn_dir, ligands_dir, status_lig: QTextEdit):
        super().__init__()
        self.codyn_dir = codyn_dir
        self.ligands_dir = ligands_dir
        self.status_lig = status_lig

    def log(self, message: str):
        # Adiciona mensagem no status_lig
        self.status_lig.append(f"[{datetime.now().strftime('%H:%M:%S')}] {message}")
        # Auto scroll
        self.status_lig.verticalScrollBar().setValue(
            self.status_lig.verticalScrollBar().maximum()
        )
    def parametrizar_ligante(self, nome_ligante: str, campo_forca: str, lig_mol2_file: str, run_folder: str):
        try:
            base_dir = os.path.join(self.codyn_dir, "BASE")
            self.log(f"Starting force field setup for {nome_ligante}, field {campo_forca}...")
            if campo_forca == "1":  # CHARMM36
                str_file = os.path.join(run_folder, f"{nome_ligante}_fix.str")
                perl_script = os.path.join(base_dir, "sort_mol2_bonds.pl")
                fix_mol2_file = os.path.join(run_folder, f"{nome_ligante}_fix.mol2")
                self.log("Correcting bonds in MOL2...")
                subprocess.run(f"perl '{perl_script}' '{lig_mol2_file}' '{fix_mol2_file}'", shell=True)
                # Avisar sobre o CGenFF online
                subprocess.run("xdg-open https://cgenff.silcsbio.com/", shell=True)
                QMessageBox.information(None, "CGenFF", f"Please generate the {nome_ligante}_fix.str file using CGenFF online before continuing.\nOpening site...")
                cgenff_script = os.path.join(base_dir, "cgenff_charmm2gmx_py3_nx2.py")
                self.log("Converting with CGenFF to GROMACS format...")
                subprocess.run(
                    f"python3 '{cgenff_script}' '{nome_ligante}' '{fix_mol2_file}' '{str_file}' charmm36-jul2021.ff",
                    shell=True,
                    cwd=run_folder
                )
                self.log("CHARMM pipeline completed.")
            elif campo_forca == "6":  # AMBER99SB
                if os.path.exists(lig_mol2_file):
                    self.log("Running ACPYPE for AMBER...")
                    subprocess.run(
                        f"acpype -i '{lig_mol2_file}' -c bcc -n 0",
                        shell=True,
                        cwd=run_folder
                    )
                    self.log("ACPYPE pipeline completed.")
                else:
                    self.log(f"Error: MOL2 file not found: {lig_mol2_file}")
            else:
                self.log(f"Force field {campo_forca} not supported.")
        except Exception as e:
            self.log(f"Error in parametrizar_ligante: {e}")

    # def convert_ligand_rdkit(self, lig_file, nome_ligante, run_folder):
    #     try:
    #         # Carregar o arquivo PDBQT como texto e extrair o bloco ATOM
    #         with open(f"{lig_file}.pdbqt", "r") as f:
    #             pdbqt_lines = f.readlines()
    #         pdb_lines = []
    #         for line in pdbqt_lines:
    #             if line.startswith("ATOM") or line.startswith("HETATM"):
    #                 # Remove as duas últimas colunas
    #                 parts = line.rstrip('\n').rsplit(None, 2)
    #                 new_line = parts[0] + '\n' if len(parts) > 0 else '\n'
    #                 pdb_lines.append(new_line)
    #         pdb_path = os.path.join(run_folder, f"{nome_ligante}.pdb")
    #         with open(pdb_path, "w") as f:
    #             f.writelines(pdb_lines)
    #         self.log("Converted PDBQT -> PDB (ATOM/HETATM block, sem carga/tipo).")

    #         # Carregar o PDB no RDKit
    #         mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
    #         if mol is None:
    #             self.log("RDKit failed to read PDB file.")
    #             return

    #         # Adicionar hidrogênios
    #         mol_h = Chem.AddHs(mol)
    #         pdb_h_path = os.path.join(run_folder, f"{nome_ligante}_h.pdb")
    #         Chem.MolToPDBFile(mol_h, pdb_h_path)
    #         self.log("Added hydrogens.")

    #         # Calcular cargas parciais de Gasteiger
    #         AllChem.ComputeGasteigerCharges(mol_h)
    #         # Salvar como MOL2
    #         mol2_path = os.path.join(run_folder, f"{nome_ligante}.mol2")
    #         Chem.MolToMol2File(mol_h, mol2_path)
    #         self.log("Generated MOL2 with charges.")

    #     except Exception as e:
    #         self.log(f"Error in RDKit ligand conversion: {e}")

    def convert_ligand_obabel(self, lig_file, nome_ligante, run_folder):
        try:
            # Converte e processa ligantes com Open Babel
            lig_file = os.path.join(self.ligands_dir, nome_ligante)
            subprocess.run(f"cp '{lig_file}.pdbqt' '{run_folder}'", shell=True)
            os.chdir(run_folder)
            self.log(f"Processing ligand {lig_file}.pdbqt...")
            subprocess.run(f"obabel -ipdbqt *.pdbqt -opdb -O {nome_ligante}.pdb", shell=True)
            self.log("Converted PDBQT -> PDB.")
            subprocess.run(f"obabel -ipdb {nome_ligante}.pdb -h -O {nome_ligante}_h.pdb", shell=True)
            self.log("Added hydrogens.")
            subprocess.run(f"obabel -ipdb {nome_ligante}_h.pdb --partialcharge gasteiger -omol2 -O {nome_ligante}.mol2", shell=True)
            self.log("Generated MOL2 with charges.")
        except Exception as e:
            self.log(f"Error in obabel ligand conversion: {e}")

    def start_preparation(self, nome_ligante, run_folder):
        try:
            self.log(f"\n######### INITIALIZING {nome_ligante} PREPARATION #########")
            # Cria pastas necessárias
            base_dir = os.path.join(self.codyn_dir, "BASE")
            runs_dir = os.path.join(self.codyn_dir, "RUNS")
            os.makedirs(self.ligands_dir, exist_ok=True)
            os.makedirs(base_dir, exist_ok=True)
            os.makedirs(runs_dir, exist_ok=True)
            self.log("Directories verified.")

            # Lê parâmetros de .form_data.txt
            data_file = os.path.join(self.codyn_dir, ".form_data.txt")
            self.log(f"Reading parameters from {data_file}...")
            with open(data_file, 'r') as f:
                lines = f.readlines()
                t_ns = int(lines[1].strip())  # Tempo em ns
                temp = int(lines[2].strip())  # Temperatura
                t_int = int(lines[10].strip())  # Integração em fs
                ff = lines[0].strip()  # Campo de força
                p_ion = lines[3].strip()  # Íon positivo
                n_ion = lines[4].strip()  # Íon negativo
            t_fs = 1_000_000 * t_ns
            t_ps = 1_000 * t_ns
            n_steps = int(t_fs / t_int)
            self.log(f"Parameters: t_ns={t_ns}, temp={temp}, t_int={t_int}, ff={ff}")
            self.log(f"Run folder is: {run_folder}")

            # Copia arquivos de apoio
            subprocess.run(f'cp -r "{os.path.join(base_dir, "charmm36-jul2021.ff")}" "{run_folder}"', shell=True)
            subprocess.run(f'cp "{os.path.join(runs_dir, f"{nome_ligante}_fix.str")}" "{run_folder}"', shell=True)
            arquivos_base = ["sort_mol2_bonds.pl", "cgenff_charmm2gmx_py3_nx2.py", "ions.mdp",
                             "em.mdp", "nvt.mdp", "npt.mdp", "md.mdp", "mmpbsa.in"]
            for fname in arquivos_base:
                src = os.path.join(base_dir, fname)
                dst = os.path.join(run_folder, fname)
                if os.path.exists(src):
                    shutil.copy2(src, dst)
                    self.log(f"Copied {fname} to run folder.")
                else:
                    self.log(f"Warning: {fname} not found in BASE.")

            lig_file = os.path.join(self.ligands_dir, nome_ligante)
            # self.convert_ligand_rdkit(lig_file, nome_ligante, run_folder)
            self.convert_ligand_obabel(lig_file, nome_ligante, run_folder)  

            # Ajustes de .mdp e .mol2
            self.log(f"Adjusting .mdp and .mol2 files for {nome_ligante}...")
            md_file = os.path.join(run_folder, "md.mdp")
            nvt_file = os.path.join(run_folder, "nvt.mdp")
            npt_file = os.path.join(run_folder, "npt.mdp")
            lig_mol2_file = os.path.join(run_folder, nome_ligante + ".mol2")
            md_nsteps = f"nsteps                  = {n_steps}  ; equivalente a {t_ps} ps ou {t_ns} ns"
            md_temp = f"ref_t                   = {temp}   {temp}                     ; reference temperature"
            lig_resid = f"{nome_ligante}"
            tc_groups = f"tc-grps                 = Protein_{nome_ligante} {p_ion}_{n_ion}_SOL       ; two coupling groups - more accurate"

            #Ajustar a linha 2 do ligante.mol2
            self.log(f"Adjusting line 2 of {lig_mol2_file} for residue ID...")
            if os.path.exists(lig_mol2_file):
                subprocess.run(f"sed -i '2s/.*/{lig_resid}/g' \"{lig_mol2_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {lig_mol2_file} não foi encontrado para ajuste da linha 2 (resid).")

            # Ajustar a linha 4 do md.mdp
            self.log(f"Adjusting line 4 of {md_file} for nsteps...")
            if os.path.exists(md_file):
                subprocess.run(f"sed -i '4s/.*/{md_nsteps}/g' \"{md_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {md_file} não foi encontrado para ajuste da linha 4 (nsteps).")

            # Ajustar a linha 32 do md.mdp
            self.log(f"Adjusting line 32 of {md_file} for tc_groups...")
            if os.path.exists(md_file):
                subprocess.run(f"sed -i '32s/.*/{tc_groups}/g' \"{md_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {md_file} não foi encontrado para ajuste da linha 32 (tc_groups).")

            # Ajustar a linha 34 do md.mdp
            self.log(f"Adjusting line 34 of {md_file} for temperature...")
            if os.path.exists(md_file):
                subprocess.run(f"sed -i '34s/.*/{md_temp}/g' \"{md_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {md_file} não foi encontrado para ajuste da linha 34 (Temperatura).")

            # Ajustar a linha 33 do nvt.mdp
            self.log(f"Adjusting line 33 of {nvt_file} for tc_groups...")
            if os.path.exists(nvt_file):
                subprocess.run(f"sed -i '33s/.*/{tc_groups}/g' \"{nvt_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {nvt_file} não foi encontrado para ajuste da linha 33 (tc_groups).")

            # Ajustar a linha 35 do nvt.mdp
            self.log(f"Adjusting line 35 of {nvt_file} for temperature...")
            if os.path.exists(nvt_file):
                subprocess.run(f"sed -i '35s/.*/{md_temp}/g' \"{nvt_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {nvt_file} não foi encontrado para ajuste da linha 35 (Temperatura).")

            # Ajustar a linha 33 do npt.mdp
            self.log(f"Adjusting line 33 of {npt_file} for tc_groups...")
            if os.path.exists(npt_file):
                subprocess.run(f"sed -i '33s/.*/{tc_groups}/g' \"{npt_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {npt_file} não foi encontrado para ajuste da linha 33 (tc_groups).")

            # Ajustar a linha 35 do npt.mdp
            self.log(f"Adjusting line 35 of {npt_file} for temperature...")
            if os.path.exists(npt_file):
                subprocess.run(f"sed -i '35s/.*/{md_temp}/g' \"{npt_file}\"", shell=True)
            else:
                QMessageBox.warning(self, "Aviso", f"O arquivo {npt_file} não foi encontrado para ajuste da linha 35 (Temperatura).")

            # Para o arquivo mol2:
            self.log(f"Adjusting ligand name in {lig_mol2_file}...")
            if os.path.exists(lig_mol2_file):
                subprocess.run(f"sed -i 's/UNL1/{nome_ligante} /g' \"{lig_mol2_file}\"", shell=True)

            self.parametrizar_ligante(nome_ligante, ff, lig_mol2_file, run_folder)

            # Gerar .gro
            gmx_path = "/usr/local/gromacs/bin/gmx"
            subprocess.run(f"{gmx_path} editconf -f {nome_ligante}_ini.pdb -o {nome_ligante}.gro", shell=True, cwd=run_folder)
            self.log(f"Ligand {nome_ligante} preparation completed successfully.")

        except Exception as e:
            self.log(f"Error during ligand prep: {e}")


