# IMPORTANDO AS BIBLIOTECAS PARA CODYN.py:
import sys
import os
import subprocess
from datetime import datetime
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QStyle,
    QTabWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QAction, QMessageBox,
    QFileDialog, QSizePolicy, QSplitter, QTextEdit, QCheckBox, QGridLayout
)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import (Qt, QTimer, QThread, pyqtSignal, QEventLoop)
from config_form import ConfiguracaoDinamica
import module_protein
import module_ligands
import module_cell
import module_equilibration_nvt
import module_equilibration_npt
import module_production

# CRIANDO A CLASSE PARA EXECUÇÃO DAS ETAPAS EM THREADS:
class RunWorkerThread(QThread):
    update_tab_signal = pyqtSignal(int)
    log_prot_signal = pyqtSignal(str)
    log_cell_signal = pyqtSignal(str)
    log_eq_signal = pyqtSignal(str)
    log_md_signal = pyqtSignal(str)
    log_lig_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)
    alert_signal = pyqtSignal(str, str)

    def __init__(self, tabs, tab_config, codyn_dir, run_folders_info, lig_dir, prot_dir):
        super().__init__()
        self.tabs = tabs
        self.tab_config = tab_config
        self.codyn_dir = codyn_dir
        self.run_folders_info = run_folders_info
        self.lig_dir = lig_dir
        self.prot_dir = prot_dir    
        # self.cancel_requested = False

    def request_cancel(self):
        self.cancel_requested = True

    def run(self):
        self.cancel_requested = False  # Reset flag at start
        self.lig_dir = self.tab_config.lig_dir.text()
        self.prot_dir = self.tab_config.prot_dir.text()
        self.update_tab_signal.emit(2) # STEP 1

        # Loop para preparo dos ligantes em cada pasta de execução:
        for info in self.run_folders_info:
            if self.cancel_requested:
                break  # Interrompe o loop se cancelado
            nome_ligante = info["nome_ligante"]
            # Passa para a classe PreparoLigantes em module_ligands.py
            prep_lig = module_ligands.PreparoLigantes(
                self.codyn_dir,
                self.lig_dir
            )
            prep_lig.log_lig_signal.connect(self.log_lig_signal.emit)  # Conecta o sinal de log
            prep_lig.alert_signal.connect(self.alert_signal.emit)  # Conecta o sinal de alerta
            prep_lig.start_preparation(nome_ligante, run_folder=info["run_folder"])

        # QTimer.singleShot(1000, lambda: self.tabs.setCurrentIndex(3))
        self.update_tab_signal.emit(3)  # STEP 2

        if self.cancel_requested:
            return  # Não continua para as próximas etapas

        # PREPARO DA PROTEÍNA:
        # Verifica extensão e prepara nome_proteina
        prot_files = [f for f in os.listdir(self.prot_dir.strip()) if os.path.isfile(os.path.join(self.prot_dir.strip(), f))]
        if len(prot_files) != 1:
            self.log_prot_signal.emit("Attention!", "It is only possible to perform the dynamics with one protein at a time.\nPlace only one file in the proteins folder.")
            self.cancel_requested = True
            return
        prot_file = prot_files[0]
        ext = os.path.splitext(prot_file)[1].lower()
        if ext not in [".pdb", ".pdbqt"]:
            self.log_prot_signal.emit("Attention!", "The protein file must have a .pdb or .pdbqt extension.")
            return
        nome_proteina = os.path.splitext(prot_file)[0]

        # Se for .pdbqt, converte para .pdb usando obabel
        prot_path = os.path.join(self.prot_dir, prot_file)
        if ext == ".pdbqt":
            pdb_path = os.path.join(self.prot_dir.strip(), f"{nome_proteina}.pdb")
            subprocess.run(f'obabel "{prot_path}" -O "{pdb_path}"', shell=True)
            os.remove(prot_path)  # Remove o arquivo .pdbqt original
            msg = f"[{datetime.now().strftime('%H:%M:%S')}] Arquivo {prot_file} convertido para {nome_proteina}.pdb usando Open Babel."
            self.log_prot_signal.emit(msg)
            prot_path = pdb_path  # Atualiza para o .pdb recém-criado

        # PREPARO DO COMPLEXO:
        # Loop para preparo do complexo em cada pasta de execução:
        for info in self.run_folders_info:
            if self.cancel_requested:
                break
            nome_proteina = info["nome_proteina"]
            nome_ligante=info["nome_ligante"]
            run_folder=info["run_folder"]
            config_data = os.path.join(self.codyn_dir, "form_data.json")
            prep_prot = module_protein.PreparoProteina(
                self.codyn_dir,
                self.prot_dir.strip(),
                self.lig_dir.strip(),
                config_data,
                run_folder,
                prot_path
            )
            prep_prot.log_prot_signal.connect(self.log_prot_signal.emit)  # Conecta o sinal de log
            prep_prot.start_preparation(nome_proteina, nome_ligante)

        if self.cancel_requested:
            return        

        # QTimer.singleShot(1000, lambda: self.tabs.setCurrentIndex(4))
        self.update_tab_signal.emit(4)  # STEP 3

        # Loop para preparo da Célula Unitária e Minimização:
        self.prep_cells = [] # Lista para armazenar instâncias de PreparoCell
        for info in self.run_folders_info:            
            if self.cancel_requested:
                break
            nome_proteina = info["nome_proteina"]
            nome_ligante = info["nome_ligante"]
            run_folder = info["run_folder"]
            config_data = os.path.join(self.codyn_dir, "form_data.json")
            prep_cell = module_cell.PreparoCell(
                self.codyn_dir,
                config_data,
                run_folder
            )
            prep_cell.log_cell_signal.connect(self.log_cell_signal.emit)
            self.prep_cells.append(prep_cell)  # Salva a referência!
            prep_cell.start_preparation(nome_proteina, nome_ligante)
            if hasattr(prep_cell, 'min_thread'):
                prep_cell.min_thread.wait()  # Espera terminar antes de seguir

        if self.cancel_requested:
            return

        # QTimer.singleShot(1000, lambda: self.tabs.setCurrentIndex(5))
        self.update_tab_signal.emit(5)  # STEP 4
        
        # Loop para equilibração NVT:
        self.prep_nvts = [] # Lista para armazenar instâncias de PrepNVT
        for info in self.run_folders_info:
            if self.cancel_requested:
                break
            nome_proteina = info["nome_proteina"]
            nome_ligante = info["nome_ligante"]
            run_folder = info["run_folder"]
            config_data = os.path.join(self.codyn_dir, "form_data.json")
            prep_nvt = module_equilibration_nvt.PrepNVT(
                self.codyn_dir,
                config_data,
                run_folder
            )
            prep_nvt.log_eq_signal.connect(self.log_eq_signal.emit)
            self.prep_nvts.append(prep_nvt)  # Salva a referência!
            prep_nvt.start_equilibration(nome_proteina, nome_ligante)
            if hasattr(prep_nvt, 'nvt_thread'):
                prep_nvt.nvt_thread.wait()  # Espera terminar antes de seguir

        if self.cancel_requested:
            return

        # Loop para equilibração NPT:
        self.prep_npts = [] # Lista para armazenar instâncias de PrepNPT
        for info in self.run_folders_info:
            if self.cancel_requested:
                break
            nome_proteina = info["nome_proteina"]
            nome_ligante = info["nome_ligante"]
            run_folder = info["run_folder"]
            config_data = os.path.join(self.codyn_dir, "form_data.json")
            prep_npt = module_equilibration_npt.PrepNPT(
                self.codyn_dir,
                config_data,
                run_folder
            )
            prep_npt.log_eq_signal.connect(self.log_eq_signal.emit)
            self.prep_npts.append(prep_npt)  # Salva a referência!
            prep_npt.start_equilibration(nome_proteina, nome_ligante)
            if hasattr(prep_npt, 'npt_thread'):
                prep_npt.npt_thread.wait()  # Espera terminar antes de seguir

        if self.cancel_requested:
            return

        # QTimer.singleShot(1000, lambda: self.tabs.setCurrentIndex(6))
        self.update_tab_signal.emit(6)  # STEP 5

        # Inicia a produção de dinâmica molecular:
        # Loop para MD Production:
        self.prep_mds = [] # Lista para armazenar instâncias de PrepMD
        for info in self.run_folders_info:
            if self.cancel_requested:
                break
            nome_proteina = info["nome_proteina"]
            nome_ligante = info["nome_ligante"]
            run_folder = info["run_folder"]
            config_data = os.path.join(self.codyn_dir, "form_data.json")
            prep_md = module_production.PrepMD(
                self.codyn_dir,
                config_data,
                run_folder
            )
            prep_md.log_md_signal.connect(self.log_md_signal.emit)
            self.prep_mds.append(prep_md)  # Salva a referência!
            
            loop = QEventLoop()
            def on_md_finished(_):
                loop.quit()
            # Conecta o sinal de término da thread de produção
            def connect_thread_and_start():
                # Garante que a thread foi criada
                if hasattr(prep_md, 'MD_Thread'):
                    prep_md.MD_Thread.finished_signal.connect(on_md_finished)
                else:
                    # Se a thread ainda não existe, conecta após a criação
                    orig_start = prep_md.start_production
                    def wrapped_start(*args, **kwargs):
                        orig_start(*args, **kwargs)
                        if hasattr(prep_md, 'MD_Thread'):
                            prep_md.MD_Thread.finished_signal.connect(on_md_finished)
                    prep_md.start_production = wrapped_start
            connect_thread_and_start()
            prep_md.start_production(nome_proteina, nome_ligante)
            loop.exec_()  # Espera até a thread terminar

        if self.cancel_requested:
            return
        
         # Muda para a aba de análise e visualização
        self.update_tab_signal.emit(7)  # STEP 6

# CRIANDO A CLASSE PRINCIPAL DA APLICAÇÃO:
class MainWindow(QMainWindow):
    # Inicializando a janela principal:
    def __init__(self):        
        super().__init__()
        self.cancel_requested = False

        # Diretórios do projeto
        self.codyn_dir = os.path.abspath(os.path.dirname(__file__))
        self.run_folder = None
        self.lig_dir = None
        self.prot_dir = None

        # Título da janela
        self.setWindowTitle("CODYN")
        
        # Tamanho da janela
        self.resize(1000, 600)
       
        # Codyn Menu Bar        
        menu_bar = self.menuBar()  
        codyn_menu = menu_bar.addMenu("Menu")        
        new_action = QAction("Configure New Run", self)
        new_action.triggered.connect(lambda: self.tabs.setCurrentIndex(1))  # Configuration tab
        codyn_menu.addAction(new_action)
        step1_action = QAction("Step 1 - Ligands Preparation", self)
        step1_action.triggered.connect(lambda: self.tabs.setCurrentIndex(2))  # Step 1
        codyn_menu.addAction(step1_action)
        step2_action = QAction("Step 2 - Complex Preparation", self)
        step2_action.triggered.connect(lambda: self.tabs.setCurrentIndex(3))  # Step 2
        codyn_menu.addAction(step2_action)
        step3_action = QAction("Step 3 - Cell Unit Preparation", self)
        step3_action.triggered.connect(lambda: self.tabs.setCurrentIndex(4))  # Step 3
        codyn_menu.addAction(step3_action)
        step4_action = QAction("Step 4 - System Equilibration", self)
        step4_action.triggered.connect(lambda: self.tabs.setCurrentIndex(5))  # Step 4
        codyn_menu.addAction(step4_action)
        step5_action = QAction("Step 5 - Molecular Dynamics Production", self)
        step5_action.triggered.connect(lambda: self.tabs.setCurrentIndex(6))  # Step 5
        codyn_menu.addAction(step5_action)
        step6_action = QAction("Step 6 - Result Analysis and Visualization", self)
        step6_action.triggered.connect(lambda: self.tabs.setCurrentIndex(7))  # Step 6
        codyn_menu.addAction(step6_action)
        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        codyn_menu.addAction(exit_action)
        
        # Help Menu        
        help_menu = menu_bar.addMenu("Help")
        about_action = QAction("About Us", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        # Tab Widget
        self.tabs = QTabWidget()
        self.tabs.tabBar().setExpanding(True)
        self.setCentralWidget(self.tabs)

        # Tabs
        self.tab_home = QWidget()
        self.tab_config = ConfiguracaoDinamica(self.codyn_dir, self.tabs)
        self.tab1 = QWidget() # Preparo Ligantes
        self.tab2 = QWidget() # Preparo do Complexo
        self.tab3 = QWidget() # Preparo Célula Unitária
        self.tab4 = QWidget() # Equilíbrio do Sistema
        self.tab5 = QWidget() # Produção
        self.tab6 = QWidget() # Análise e Visualização

        # Add tabs to the tab widget
        self.tabs.addTab(self.tab_home, "HOME")
        self.tabs.addTab(self.tab_config, "CONFIG")
        self.tabs.addTab(self.tab1, "STEP 1")
        self.tabs.addTab(self.tab2, "STEP 2")
        self.tabs.addTab(self.tab3, "STEP 3")
        self.tabs.addTab(self.tab4, "STEP 4")
        self.tabs.addTab(self.tab5, "STEP 5")
        self.tabs.addTab(self.tab6, "STEP 6")

        # Initialize tabs
        self.init_tab_home()
        self.tab_config.saved.connect(self.start_config_run_thread)
        self.init_tab1()
        self.init_tab2()
        self.init_tab3()
        self.init_tab4()
        self.init_tab5()
        self.init_tab6()

    def init_tab_home(self):
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignCenter)
        
        layout.addStretch()  # Espaço flexível no topo para centralização vertical
        
        # Logo
        logo_label = QLabel()
        logo_path = os.path.join(self.codyn_dir, "ICONS", "LOGO_CODYN2.png")
        if os.path.exists(logo_path):
            pix = QPixmap(logo_path).scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            logo_label.setPixmap(pix)
        else:
            logo_label.setText("Logo not found")
        logo_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(logo_label)        
        layout.addSpacing(10)  # Pequeno espaço fixo entre logo e título
        
        # Título
        title = QLabel("Automated Molecular Dynamics")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size:18pt;font-weight:bold;")
        layout.addWidget(title)

        # Botão Start abaixo do splitter
        layout.addSpacing(10)
        start_btn = QPushButton("START")
        start_btn.setFixedWidth(150)
        start_btn.setStyleSheet("font-size:12pt;font-weight:bold;")
        start_btn.clicked.connect(lambda: self.tabs.setCurrentIndex(1))
        layout.addWidget(start_btn, alignment=Qt.AlignCenter)

        # Botão Exit centralizado logo abaixo do botão Start
        layout.addSpacing(10)
        exit_btn = QPushButton("EXIT")
        exit_btn.setFixedWidth(150)
        exit_btn.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        exit_btn.clicked.connect(self.close)
        layout.addWidget(exit_btn, alignment=Qt.AlignCenter)

        layout.addStretch()  # Espaço flexível em baixo para centralização vertical

        self.tab_home.setLayout(layout)

    def init_tab1(self):
        main_layout = QVBoxLayout()
        label = QLabel("Ligands Preparation")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("font-size:14pt;font-weight:bold;")
        main_layout.addWidget(label)
        main_layout.addSpacing(10)  # Espaço fixo entre o título e o status

        # Status do preparo dos ligantes
        status_lig = QTextEdit()
        status_lig.setReadOnly(True)
        status_lig.setStyleSheet("background-color: #000000; font-family: monospace; color: #FFFFFF;")
        self.status_lig = status_lig
        main_layout.addWidget(status_lig)

        # Botão Cancelar
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(self.cancel_all_processes)
        main_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)

        self.tab1.setLayout(main_layout)

    def init_tab2(self):
        main_layout = QVBoxLayout()
        label = QLabel("Complex Preparation")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("font-size:14pt;font-weight:bold;")
        main_layout.addWidget(label)
        main_layout.addSpacing(10)  # Espaço fixo entre o título e o status

        # Status do preparo do complexo
        status_prot = QTextEdit()
        status_prot.setReadOnly(True)
        status_prot.setStyleSheet("background-color: #000000; font-family: monospace; color: #FFFFFF;")
        self.status_prot = status_prot
        main_layout.addWidget(status_prot)

        # Botão Cancelar
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(self.cancel_all_processes)
        main_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)

        self.tab2.setLayout(main_layout)
    
    def init_tab3(self):
        main_layout = QVBoxLayout()
        label = QLabel("Unit Cell Preparation and Minimization")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("font-size:14pt;font-weight:bold;")
        main_layout.addWidget(label)
        main_layout.addSpacing(10)  # Espaço fixo entre o título e o status

        # Status do preparo da célula unitária
        status_cell = QTextEdit()
        status_cell.setReadOnly(True)
        status_cell.setStyleSheet("background-color: #000000; font-family: monospace; color: #FFFFFF;")
        self.status_cell = status_cell
        main_layout.addWidget(status_cell)

        # Botão Cancelar
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(self.cancel_all_processes)
        main_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)

        self.tab3.setLayout(main_layout)

    def init_tab4(self):
        main_layout = QVBoxLayout()
        label = QLabel("System Equilibration")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("font-size:14pt;font-weight:bold;")
        main_layout.addWidget(label)
        main_layout.addSpacing(10)  # Espaço fixo entre o título e o status

        # Status do equilíbrio do sistema
        status_eq = QTextEdit()
        status_eq.setReadOnly(True)
        status_eq.setStyleSheet("background-color: #000000; font-family: monospace; color: #FFFFFF;")
        self.status_eq = status_eq
        main_layout.addWidget(status_eq)

        # Botão Cancelar
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(self.cancel_all_processes)
        main_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)

        self.tab4.setLayout(main_layout)

    def init_tab5(self):
        main_layout = QVBoxLayout()
        label = QLabel("Molecular Dynamics Production")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("font-size:14pt;font-weight:bold;")
        main_layout.addWidget(label)
        main_layout.addSpacing(10)  # Espaço fixo entre o título e o status

        # Status da produção
        status_prod = QTextEdit()
        status_prod.setReadOnly(True)
        status_prod.setStyleSheet("background-color: #000000; font-family: monospace; color: #FFFFFF;")
        self.status_prod = status_prod
        main_layout.addWidget(status_prod)

        # Botão Cancelar
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(self.cancel_all_processes)
        main_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)

        self.tab5.setLayout(main_layout)
    def init_tab6(self):        
        main_layout = QVBoxLayout()

        # Label centralizada no topo
        label = QLabel("Select the results you want to view:")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("font-size:14pt;font-weight:bold;")
        main_layout.addWidget(label)
        main_layout.addSpacing(1)

        # Grid para duas colunas de opções
        grid = QGridLayout()
        options = [ "PotentialEnergy-EM", "Temperature-NVT", "Pressure-NPT", "Density-NPT", "TotalEnergy-NPT", "PROT-RMSD", "LIG-RMSD", "PROT_LIG-RMSD", "PROT-RMSF", "PROT_LIG-HBOND", "PROT_LIG-COM", "PROT_LIG-SASA"]
        self.checkboxes_tab6 = {}

        # Distribui as opções em duas colunas
        for i, opt in enumerate(options):
            checkbox = QCheckBox(opt)
            self.checkboxes_tab6[opt] = checkbox
            row = i % 2
            col = i // 2
            grid.addWidget(checkbox, row, col)

        # Centraliza o grid
        grid_widget = QWidget()
        grid_widget.setLayout(grid)
        grid_layout = QHBoxLayout()

        grid_layout.addWidget(grid_widget)
        # grid_layout.addStretch()
        main_layout.addLayout(grid_layout)
        # main_layout.addSpacing(20)

        # Janela para exibição de gráficos (placeholder)
        self.graph_area = QTextEdit()                      
        self.graph_area.setReadOnly(True)
        self.graph_area.setFixedHeight(350)
        self.graph_area.setStyleSheet("background-color: #222; color: #fff; font-family: monospace;")
        main_layout.addWidget(self.graph_area)

        # Botões
        # Layout horizontal para os botões
        button_layout = QHBoxLayout()

        # Botão RUN centralizado
        run_btn = QPushButton("RUN")
        run_btn.setFixedWidth(150)
        run_btn.setStyleSheet("font-size:12pt;font-weight:bold;background-color:green;color:white;")
        # run_btn.clicked.connect(self.sua_funcao_de_analise)  # Conecte à função desejada
        button_layout.addWidget(run_btn, alignment=Qt.AlignCenter)

        # Botão Cancelar
        btn_cancel = QPushButton("CANCEL")
        btn_cancel.setFixedWidth(150)
        btn_cancel.setStyleSheet("font-size:12pt;font-weight:bold;background-color:red;color:white;")
        btn_cancel.clicked.connect(self.cancel_all_processes)
        button_layout.addWidget(btn_cancel, alignment=Qt.AlignCenter)

        main_layout.addLayout(button_layout)
        self.tab6.setLayout(main_layout)

    # Método para cancelar todos os processos e voltar para a aba HOME
    def cancel_all_processes(self):
        reply = QMessageBox.question(
            self,
            "Attention",
            "Do you really want to cancel all processes\n and return to the HOME tab?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            if hasattr(self, 'worker'):
                self.worker.request_cancel()
            self.tabs.setCurrentIndex(0)

    def show_alert(self, title, message):
        QMessageBox.information(self, title, message)

    def log_status_lig(self, message):
        self.status_lig.append(message)
        self.status_lig.verticalScrollBar().setValue(self.status_lig.verticalScrollBar().maximum())

    def log_status_prot(self, message):
        # Exibe no QTextEdit ou, se necessário, chama QMessageBox
        if message.startswith("Attention!"):
            QMessageBox.warning(self, "Attention!", message[9:])
        else:
            self.status_prot.append(message)
            self.status_prot.verticalScrollBar().setValue(self.status_prot.verticalScrollBar().maximum())

    def log_status_cell(self, message):
        self.status_cell.append(message)
        self.status_cell.verticalScrollBar().setValue(self.status_cell.verticalScrollBar().maximum())

    def log_status_eq(self, message):
        self.status_eq.append(message)
        self.status_eq.verticalScrollBar().setValue(self.status_eq.verticalScrollBar().maximum())

    def log_status_md(self, message):
        self.status_prod.append(message)
        self.status_prod.verticalScrollBar().setValue(self.status_prod.verticalScrollBar().maximum())

    def start_config_run_thread(self, run_folders_info):
        self.worker = RunWorkerThread(
            self.tabs,
            self.tab_config,
            self.codyn_dir,
            run_folders_info,
            self.tab_config.lig_dir.text(),
            self.tab_config.prot_dir.text()
        )
        self.worker.alert_signal.connect(self.show_alert, Qt.BlockingQueuedConnection)
        self.worker.log_prot_signal.connect(self.log_status_prot)
        self.worker.log_cell_signal.connect(self.log_status_cell)
        self.worker.log_eq_signal.connect(self.log_status_eq)
        self.worker.log_md_signal.connect(self.log_status_md)
        self.worker.log_lig_signal.connect(self.log_status_lig)
        self.worker.update_tab_signal.connect(self.tabs.setCurrentIndex)
        self.worker.finished_signal.connect(lambda _: QMessageBox.information(self, "Concluído", "Dinâmica finalizada!"))
        self.worker.start()

    def closeEvent(self, event):
        if hasattr(self, 'worker') and self.worker.isRunning():
            self.worker.request_cancel()
            self.worker.wait()
        event.accept()

    def show_about(self):
        QMessageBox.about(
            self,
            "CODYN",
            "                          An Open Source\n   Automated Molecular Dynamics Tool\n   \nVersion 1.0\n   \n© July-2025\n   \nDeveloped by:\n   Moisés Maia Neto and Gustavo Scheiffer\n   Federal University of Paraná (UFPR) - Brazil\n \nContact:\n   moimaian@gmail.com"
        )


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
