import os
from PyQt5.QtWidgets import (
    QWidget, QLabel, QFileDialog, QLineEdit, QPushButton, QFormLayout, QVBoxLayout, QHBoxLayout, QMessageBox, QComboBox
)
from PyQt5.QtCore import Qt, pyqtSignal
from datetime import datetime
import json

class ConfiguracaoDinamica(QWidget):
    saved = pyqtSignal(list)         # sinal que será emitido ao salvar

    def __init__(self, codyn_dir, tabs):
        super().__init__()
        self.codyn_dir = codyn_dir
        self.tabs = tabs
        self.setWindowTitle("CODYN - Molecular Dynamics Configuration")
        self.setGeometry(200, 200, 500, 600)

        # Layout principal
        layout = QVBoxLayout()
        form_layout = QFormLayout()

        # -------- Campo de Força como ComboBox --------
        self.ff = QComboBox()
        self.ff.addItems([
            "1 - Charmm36_2021",
            "6 - AMBER99SB",
            "9 - Charmm27",
            "15 - GROMOS96_54a7",
            "16 - OPLS_AA/L"
        ])
        form_layout.addRow("Force Field:", self.ff)

        # -------- Campos de texto normais --------
        self.t_ns_md = QLineEdit("100")
        self.t_ns_eq = QLineEdit("5")
        self.temp = QLineEdit("310")
        self.p_ion = QLineEdit("SOD")
        self.n_ion = QLineEdit("CLA")
        self.conc_ions = QLineEdit("0.15")
        self.t_int = QLineEdit("2")        
        self.n_threads = QLineEdit(str(os.cpu_count()))

        # -------- Modelo de Solvente como ComboBox --------
        self.mod_sol = QComboBox()
        self.mod_sol.addItems([
            "1 - TIP3P",
            "2 - TIP4P",
            "3 - TIP5P",
            "4 - SPC",
            "5 - SPC/E"
        ])

        # -------- Adicionando campos ao formulário --------
        form_layout.addRow("Production time (ns):", self.t_ns_md)
        form_layout.addRow("Equilibration time (ns):", self.t_ns_eq)
        form_layout.addRow("Integration Time (fs):", self.t_int)
        form_layout.addRow("Temperature (K):", self.temp)
        form_layout.addRow("Positive Ion (M+):", self.p_ion)
        form_layout.addRow("Negative Ion (A-):", self.n_ion)
        form_layout.addRow("Ion Concentration (M):", self.conc_ions)
        form_layout.addRow("Solvent Model:", self.mod_sol)
        form_layout.addRow("Number of Threads:", self.n_threads)

        # -------- Diretórios de Ligantes e Proteínas --------
        self.lig_dir = QLineEdit()
        self.lig_dir.setReadOnly(True)
        self.lig_dir.setText(os.path.join(self.codyn_dir, "LIGANDS"))
        lig_btn = QPushButton("Search")
        lig_btn.clicked.connect(self.select_lig_dir)
        lig_layout = QHBoxLayout()
        lig_layout.addWidget(self.lig_dir)
        lig_layout.addWidget(lig_btn)

        self.prot_dir = QLineEdit()
        self.prot_dir.setReadOnly(True)
        self.prot_dir.setText(os.path.join(self.codyn_dir, "TARGETS"))
        prot_btn = QPushButton("Search")
        prot_btn.clicked.connect(self.select_prot_dir)
        prot_layout = QHBoxLayout()
        prot_layout.addWidget(self.prot_dir)
        prot_layout.addWidget(prot_btn)

        form_layout.addRow("Ligands Folder:", lig_layout)
        form_layout.addRow("Target Folder:", prot_layout)

        # Adiciona o layout do formulário ao layout principal
        layout.addLayout(form_layout)

        # Botões
        # Layout horizontal para os botões
        button_layout = QHBoxLayout()

        # Botão Cancelar - à esquerda
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(lambda: self.tabs.setCurrentIndex(0))  # Volta para a aba HOME
        button_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)        

        # Botão Salvar:
        btn_save = QPushButton("RUN")
        btn_save.setFixedWidth(150)
        btn_save.setStyleSheet("font-size:12pt;font-weight:bold;background-color:green;color:white;")
        self.btn_save = btn_save
        self.btn_save.clicked.connect(self.salvar_dados)
        button_layout.addWidget(btn_save, alignment=Qt.AlignCenter)

        # Adiciona o layout de botões ao layout principal
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def select_lig_dir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Ligands Folder")
        if directory:
            self.lig_dir.setText(directory)

    def select_prot_dir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Target Folder")
        if directory:
            self.prot_dir.setText(directory)

    def salvar_dados(self):
        try:
            # Pegando apenas o número da escolha nos ComboBox
            ff_value = self.ff.currentText().split('-')[0].strip()
            mod_sol_value = self.mod_sol.currentText().split('-')[0].strip()
            lig_dir = self.lig_dir.text().strip()
            prot_dir = self.prot_dir.text().strip()

            form_data = {
                "force_field": ff_value,
                "production_time_ns": self.t_ns_md.text(),
                "equilibration_time_ns": self.t_ns_eq.text(),
                "integration_time_fs": self.t_int.text(),
                "temperature_K": self.temp.text(),
                "positive_ion": self.p_ion.text(),
                "negative_ion": self.n_ion.text(),
                "ion_concentration_M": self.conc_ions.text(),
                "solvent_model": mod_sol_value,
                "n_threads": self.n_threads.text(),
                "ligands_dir": lig_dir,
                "protein_dir": prot_dir
            }
            # Verifica se todos os campos estão preenchidos
            if any(not item for item in form_data):
                QMessageBox.warning(self, "Attention!", "There are empty configuration fields.")
                self.tabs.setCurrentIndex(1)  # Reinicia a aba de configuração
                return

            # Salvar em arquivo oculto:
            file_path = os.path.join(self.codyn_dir, "form_data.json")
            with open(file_path, "w") as f:
                json.dump(form_data, f, indent=4)
            
            # Emite o log de salvamento:
            QMessageBox.information(self, "Sucess!", f"Parameters saved in: {file_path}")

            # Lista ligantes e proteínas:            
            lig_files = [f for f in os.listdir(lig_dir) if os.path.isfile(os.path.join(lig_dir, f))]
            prot_files = [f for f in os.listdir(prot_dir) if os.path.isfile(os.path.join(prot_dir, f))]

            # Verifica se há arquivos de ligantes
            if not lig_files:
                QMessageBox.warning(self, "Attention", "No ligands found in the specified folder.")
                self.tabs.setCurrentIndex(0)
                return
            
            # Verifica se há apenas um arquivo de proteína
            if len(prot_files) != 1:
                QMessageBox.warning(self, "Atenção", "Só é possível realizar a dinâmica com uma proteína por vez.\nColoque apenas um arquivo na pasta de proteínas.")
                self.tabs.setCurrentIndex(0)
                return

            # Cria pastas de execução:
            data_str = datetime.now().strftime("%Y-%m-%d_%H-%M")
            runs_dir = os.path.join(self.codyn_dir, "RUNS")
            os.makedirs(runs_dir, exist_ok=True)
            run_folders_info = []
            for prot_file in prot_files:
                nome_proteina = os.path.splitext(prot_file)[0]
                for lig_file in lig_files:
                    nome_ligante = os.path.splitext(lig_file)[0]
                    run_folder = os.path.join(runs_dir, f"{data_str}_{nome_proteina}_{nome_ligante}")
                    os.makedirs(run_folder, exist_ok=True)
                    run_folders_info.append({
                        "run_folder": run_folder,
                        "nome_proteina": nome_proteina,
                        "nome_ligante": nome_ligante
                    })
            self.run_folders_info = run_folders_info

            QMessageBox.information(self, "Sucesso", "Pastas criadas:\n" + "\n".join([rf["run_folder"] for rf in run_folders_info]))
            self.saved.emit(run_folders_info)

        except Exception as e:
            QMessageBox.critical(self, "Erro", f"Erro ao salvar o arquivo: {str(e)}")
