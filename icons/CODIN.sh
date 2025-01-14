#!/bin/bash
#################################################################################################################################
#                                               CODIN VERSION 1.1.2:                                                             #
#                                      Developed by Moisés Maia Neto - 07/2024                                                  #
#################################################################################################################################
# Identifying the directory where the CODIN is located and checking the TARGETS and LIGANDS directories.
CURRENT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd $CURRENT_DIR
CODIN_DIR=$(pwd)
if [ -d "LIGANDS" ]; then
    echo "LIGANDS FOLDER IS PRESENT!"
    else
        dialog --title "ATTENTION" --msgbox "THE DIRECTORY CONTAINING THE TARGETS IS MISSING. CREATE A DIRECTORY CALLED "LIGANDS"!" 10 60 ;
    fi
if [ -d "TARGETS" ]; then
    echo "TARGETS FOLDER IS PRESENT!"
    else
        dialog --title "ATTENTION" --msgbox "THE DIRECTORY CONTAINING THE LIGANDS IS MISSING. CREATE A DIRECTORY CALLED "TARGETS"!" 10 60 ;
fi

#Atualizando o executável Gromacs:
source /usr/local/gromacs/bin/GMXRC

#Variáveis de acesso global predefinidas:
base="$CODIN_DIR/BASE"
ligands="$CODIN_DIR/LIGANDS"
din="$CODIN_DIR"
ff="1"
t_ns="100"
temp="310"
p_ion="SOD"
n_ion="CLA"
conc_ions="0.15"
mod_sol="1"
box_x="10.00000"
box_y="10.00000"
box_z="10.00000"
t_int="2"
n_threads="16"

#Função de acesso global para remoção de parâmetros:
remover_parametros(){
rm /$CODIN_DIR/.form_data.txt
rm /$CODIN_DIR/.menu_escolha.txt
}

# Função para exibir o menu principal
programa=("- C O D I N - Dinâmica Multitarget com Gromacs 2023.4")
show_main_menu() {
    dialog --clear --backtitle "$programa" \
    --title "MENU PRINCIPAL" \
    --menu "Selecione uma opção:" 15 60 5 \
    1 "Configurações da Dinâmica" \
    2 "Preparo dos ligantes" \
    3 "Equilibração(NVT/NPT) & Produção(MD)" \
    4 "Energia de Ligação(MMGBSA)" \
    5 "SAIR" \
    2> .menu_escolha.txt

    opcao=$(cat .menu_escolha.txt)
    case $opcao in
        1) show_form ;;
        2) executa_preparo_ligantes ; dialog --msgbox "Verifique a pasta $ligands." 0 0 ;;
        3) executa_producao ;;
        4) executa_mmgbsa ;;
        5) remover_parametros; clear ; exit 0 ;;
        *) remover_parametros; clear ; exit 0 ;;
    esac
}

###########################################################
# Função para configurar parâmetros e variáveis da dinâmica:
###########################################################
show_form() {
    form_data
    form_data=$(dialog --backtitle "$programa" \
    --title "OPÇÃO 1: PARÂMETROS DE CONFIGURAÇÃO" \
    --form "Preencha os campos abaixo:" 35 95 0 \
    "Campo de Força:" 2 4 "$ff" 2 35 20 0 \
    "1-Charmm36_2021 / 6-AMBER99SB / 9-Charmm27 / 15-GROMOS96_54a7 / 16-OPLS_AA/L" 3 4 "" 0 0 0 0 \
    "Tempo da simulação (ns):" 5 4 "$t_ns" 5 35 20 0 \
    "Temperatura do sistema (K):" 7 4 "$temp" 7 35 20 0 \
    "Íon Positivo:" 9 4 "$p_ion" 9 35 20 0 \
    "SOD-Sódio / POT-Potássio / MG-Magnésio / CAL-Cálcio / FE2P-Ferroso / FE3P-Férrico" 10 4 "" 0 0 0 0 \
    "Íon Negativo:" 12 4 "$n_ion" 12 35 20 0 \
    "CLA-Cloreto / OH-Hidróxido" 13 4 "" 0 0 0 0 \
    "Concentração dos Íons (M):" 15 4 "$conc_ions" 15 35 20 0 \
    "Modelos de solvente:" 17 4 "$mod_sol" 17 35 20 0 \
    "1-TIP3P / 2-TIP4P / 3-TIP5P / 4-SPC / 5-SPC/E" 18 4 "" 0 0 0 0 \
    "Box_X_Size (nm):" 20 4 "$box_x" 20 35 20 0 \
    "Box_Y_Size (nm):" 22 4 "$box_y" 22 35 20 0 \
    "Box_Z_Size (nm):" 24 4 "$box_z" 24 35 20 0 \
    "Tempo de Integração (fs):" 26 4 "$t_int" 26 35 20 0 \
    "Número de Threads (CPU):" 28 4 "$n_threads" 28 35 20 0 \
    3>&1 1>&2 2>&3)

    # Verificar se o usuário cancelou o formulário
    if [ $? -eq 0 ]; then
        # Atribuir os dados inseridos às variáveis
        echo "$form_data" > .form_data.txt
        source .form_data.txt
        #show_main_menu
    else
        echo "Cancelado."
        show_main_menu
    fi
}
show_main_menu

###########################################################
# Função para configurar parâmetros e variáveis da dinâmica
###########################################################
executa_preparo_ligantes() {
#Atualiza as variáveis:
data=/$CODIN_DIR/.form_data.txt
ff=$(sed -n '1p' $data)
t_ns=$(sed -n '2p' $data)
temp=$(sed -n '3p' $data)
p_ion=$(sed -n '4p' $data)
n_ion=$(sed -n '5p' $data)
conc_ions=$(sed -n '6p' $data)
mod_sol=$(sed -n '7p' $data)
box_x=$(sed -n '8p' $data)
box_y=$(sed -n '9p' $data)
box_z=$(sed -n '10p' $data)
t_int=$(sed -n '11p' $data)
n_threads=$(sed -n '12p' $data)

echo "CÁLCULO DE DINÂMICA MOLECULAR PARA MÚLTIPLOS LIGANTES"
echo "GROMACS 2023.4"
echo " "
echo "PARAMETRIZAÇÃO DOS LIGANTES"
echo " "
date +"%d-%m-%Y %H:%M:%S"
sleep 1

#Declarando endereços e outras variáveis úteis para configurações da dinâmica:
base="$CODIN_DIR/BASE"
ligands="$CODIN_DIR/LIGANDS"
din="$CODIN_DIR"
t_fs=$((1000000 * t_ns)) # Tempo de dinâmica em fentossegundos (10e-15 s)
t_ps=$((1000 * t_ns)) # Tempo de dinâmica em Picossegundos (10e-12 s)
n_steps=$((t_fs / "$t_int")) #Número de etapas(nsteps) que deve ser incluído no arquivo md.mdp 

#PREPARO DO LIGANTE:
echo "####### PREPARO DE LIGANTES #######"
echo " "
echo "Transfira os arquivos .pdbqt com as melhores poses dos ligantes..."
echo "...para a pasta $ligands"
echo " "
echo "ATENÇÃO! RENOMEIE COM UMA SIGLA DE 3 LETRAS PARA QUE ELA SEJA RECONHECIDA COMO RESÍDUO"
echo " "
echo "Quando terminar digite qualquer tecla para continuar"
echo " "
read -n 1

# Entre na pasta:
echo "Entrando na pasta de ligantes"
echo " "
cd $ligands

# Converter para .pdb:
echo "Convertendo arquivos .pdbqt para .pdb"
echo " "
sleep 1
obabel -ipdbqt *.pdbqt -opdb -O *.pdb 

# Adicionar os hidrogênios:
echo " "  
echo "Adicionando todos os hidrogênios faltantes"
echo " "
sleep 1
obabel -ipdb *.pdb -h -O*_h.pdb 

# Adicionar cargas parciais e converter para mol2:
echo " "  
echo "Adicionando cargas parciais por gasteiger e convertendo para sybil .mol2"
echo " "
sleep 1
obabel -ipdb *_h.pdb --partialcharge gasteiger -omol2 -O *.mol2

# Remover os demais arquivos de ligantes que não serão usados:
echo " "  
echo "Removendo arquivos desnecessários!..."
echo " "
rm *.pdbqt
rm *.pdb
sleep 1

for lig in *.mol2; do
    L=$(basename "$lig" _h.mol2)
    l="${L,,}"
    mkdir -p "$L"
    cp "$lig" "$ligands"/"$L"
    cp "$base"/sort_mol2_bonds.pl "$ligands"/"$L"
    cp "$base"/cgenff_charmm2gmx_py3_nx2.py "$ligands"/"$L"
    cp "$base"/ions.mdp "$ligands"/"$L"
    cp "$base"/em.mdp "$ligands"/"$L"
    cp "$base"/nvt.mdp "$ligands"/"$L"
    cp "$base"/npt.mdp "$ligands"/"$L"
    cp "$base"/md.mdp "$ligands"/"$L"
    cp "$base"/mmpbsa.in "$ligands"/"$L"
    cp -r $base/charmm36-jul2021.ff "$ligands"/"$L"
    
    rm $lig
    cd $ligands/$L
    #AJUSTANDO AS CONFIGURAÇÕES DOS ARQUIVOS .mdp
    sed -i "s/Protein_LIG/Protein_"$L"/g" nvt.mdp
    sed -i "s/Protein_LIG/Protein_"$L"/g" npt.mdp
    sed -i "s/Protein_LIG/Protein_"$L"/g" md.mdp  
    md_nsteps="nsteps                  = $n_steps  ; equivalente a $t_ps ps ou $t_ns ns"
    md_temp="ref_t                   = $temp   $temp                     ; reference temperature, one for each group, in K"
    sed -i "4s/.*/$md_nsteps/g" md.mdp
    sed -i "34s/.*/$md_temp/g" md.mdp
    sed -i "35s/.*/$md_temp/g" nvt.mdp
    sed -i "35s/.*/$md_temp/g" npt.mdp
    
    # Corrigindo as ligações:
    echo " "  
    echo "Corrigindo a ordem de ligações dos ligantes..."
    echo " "
    sed -i "s/UNL1/$L /g" "$lig"
    sed -i "s/"$L"_h\.pdb/$L/g" "$lig"
    perl sort_mol2_bonds.pl "$lig" "$L"_fix.mol2
    
    # Preparando o arquivo .srt:
    echo " "  
    echo "Prepare o arquivo "$L"_fix.str do ligante a partir do arquivo "$L"_fix.mol2 usando o Cgenff"
    echo " "
    echo "Abrindo o portal Cgenff..."
    echo " "
    echo "Quando terminar digite qualquer tecla para continuar"
    echo " "  
    sleep 1
    xdg-open https://cgenff.silcsbio.com/
   
    read -n 1
    
    # Preparando o arquivo _ini.pdb:
    echo " "    
    echo "Preparando o arquivo "$l"_ini.pdb do ligante a partir do arquivo "$L"_fix.str"
    echo "ATENÇÃO! O SCRIPT EM PYTHON GERA EM MINÚSCULO O NOME "$l"_ini.pdb do ligante..."
    echo "...ISSO INFLUENCIA NO RECONHECIMENTO DO COMANDO GMX EDITCONF PARA GERAR $l.gro"
    echo " "
    sleep 1
    
    python3 cgenff_charmm2gmx_py3_nx2.py "$L" "$L"_fix.mol2 "$L"_fix.str charmm36-jul2021.ff

    # Gerando o arquivo .gro:
    echo " "    
    echo "Gerando o arquivo "$l".gro do ligante"
    echo " "
    sleep 1    
    gmx editconf -f "$l"_ini.pdb -o "$l".gro
    cd ../
done
cd $din
}
show_main_menu

##############################################################
# FUNÇÃO PARA EXECUTAR A EQUILIBRAÇÃO E A PRODUÇÃO DA DINÂMICA:
##############################################################
executa_producao() {
# Início do contador de tempo
start_time=$(date +%s)
source /usr/local/gromacs/bin/GMXRC

#Atualiza as variáveis:
base="$CODIN_DIR/BASE"
ligands="$CODIN_DIR/LIGANDS"
din="$CODIN_DIR"
targets="$CODIN_DIR/TARGETS"
data=/$CODIN_DIR/.form_data.txt
ff=$(sed -n '1p' $data)
t_ns=$(sed -n '2p' $data)
temp=$(sed -n '3p' $data)
p_ion=$(sed -n '4p' $data)
n_ion=$(sed -n '5p' $data)
conc_ions=$(sed -n '6p' $data)
mod_sol=$(sed -n '7p' $data)
box_x=$(sed -n '8p' $data)
box_y=$(sed -n '9p' $data)
box_z=$(sed -n '10p' $data)
t_int=$(sed -n '11p' $data)
n_threads=$(sed -n '12p' $data)

echo " "
echo "####### PREPARO DA PROTEINA #######"
echo " "
echo "Transfira o arquivo proteina.pdb para a pasta $din/TARGETS"
echo "ATENÇÃO AO HIDROGÊNIOS TERMINAIS. PODE SER NECESSÁRIO USAR O CHARMM-GUI PARA CORREÇÕES"
echo " "
echo "Quando terminar digite qualquer tecla para continuar"
echo " "
read -n 1

cd $ligands

# Copiando o arquivo da proteína .pdb para as pastas dos ligantes:
echo " "
echo "Copiando o arquivo da proteína para as pastas dos ligantes"
echo " "
for lig in "$ligands"/*/; do
    L=$(basename "$lig")
    l="${L,,}"
    cp "$targets/proteina.pdb" "$ligands/$L/"
    cd "$ligands/$L"

    #Prepare a topologia da proteína:
    echo " "
    echo "Preparando a Topologia da proteína..."
    echo "Selecionando campo de força CHARMM36_JUL2021.FF e Modelo de Água TIP3P"
    echo " "
    gmx pdb2gmx -f proteina.pdb -o proteina.gro<<EOF
1
$mod_sol
EOF
    
    # Atualize o arquivo topol.top:
    echo "Referenciando os arquivos de parametrização $l.prm e de topologia $l.itp no arquivo topol.top:"
    echo ""

    t1_1="; Include ligand parameters"
    t1_2="#include \""$l".prm\""

    t2_1="; Include ligand topology"
    t2_2="#include \""$l".itp\""

    t3="$L     1"
    sed -i -e "/\[ moleculetype \]/i$t1_1\n$t1_2\\"$'\n' -e "/; Include water topology/i$t2_1\n$t2_2\\"$'\n' topol.top
    echo "$t3" >> topol.top
    sed -i 's/Protein_chain_A/Protein        /g' topol.top

    # PREPARO DO ARQUIVO complexo.gro:
    echo ""
    echo "Prepare o arquivo complexo.gro adicionando o "$l".gro em proteina.gro"
    echo "Atualize o número total de átomos e ajuste a posição do LIG na sequencia"
    echo ""
    # Declare a variável box que define o tamanho da caixa de dinâmica:
    box="   "$box_x"   "$box_y"  "$box_z" "    
    # Remova as linhas 1, 2 e a última do arquivo $l.gro e salve como $l_t.gro:
    sed '1d;2d;$d' $l.gro > "$l"_t.gro
    # Remova a última linha do arquivo proteina.gro e salve como complexo.gro:
    sed '$d' proteina.gro > complexo.gro
    # Copie o conteúdo de $l_t.gro para complexo.gro:
    cat "$l"_t.gro >> complexo.gro
    # Extrair valores numéricos da linha 2 dos arquivos
    atomos_ligante=$(awk 'NR==2 {print $1}' "$l.gro")
    atomos_proteina=$(awk 'NR==2 {print $1}' proteina.gro)
    # Somar os valores
    atomos_complexo=$(awk "BEGIN {printf \"%.0f\", $atomos_ligante + $atomos_proteina}")
    # Substituir a linha 2 do arquivo complexo.gro com o novo valor numérico
    awk -v atomos_complexo="$atomos_complexo" 'NR==2 {$1=atomos_complexo}1' complexo.gro > complexo_temp.gro
    mv complexo_temp.gro complexo.gro
    rm complexo_temp.gro
    # Copie o tamanho da caixa $box para complexo.gro:
    echo "$box" >> complexo.gro

    #PREPARO A CÉLULA DE SOLVATAÇÃO:
    #CONFIRA AS CONFIGURAÇÕES DO ARQUIVO ions.mdp
    echo ""
    echo "Preparando a célula unitária de solvatação..."
    echo " "
    sleep 1
    gmx editconf -f complexo.gro -o newbox.gro -bt cubic -c -d 1.0
    wait
    gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

    echo ""
    echo "Adicionando íons e ajustando o balanço de cargas..."
    echo " "
    sleep 1
    gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
    wait
    gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname "$p_ion" -nname "$n_ion" -neutral -conc "$conc_ions" <<EOF
    15
EOF
    echo "CASO APRESENTE ERRO ENTRE NO ARQUIVO IONS.ITP, NA PASTA DO CAMPO DE FORÇA, E VERIFIQUE SE A SIGLA DO ÍON ESTÁ CORRETA."
    
    #REALIZANDO A MINIMIZAÇÃO DE ENERGIA:
    #CONFIRA AS CONFIGURAÇÕES DO ARQUIVO em.mdp
    echo "Realizando a Minimização de Energia restrita..."
    echo " "
    sleep 1
    gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
    wait
    gmx mdrun -v -deffnm em

    echo "Gerando o gráfico de Energia Potencial..."
    echo " "
    sleep 1
    gmx energy -f em.edr -o potential.xvg <<EOF
    11
    0
EOF

    echo "Gerando o grupo de índices..."
    echo " "
    sleep 1
    gmx make_ndx -f "$l".gro -o index_"$l".ndx  <<EOF
    0 & ! a H*
    q
EOF

    gmx genrestr -f "$l".gro -n index_"$l".ndx -o posre_"$l".itp -fc 1000 1000 1000 <<EOF
    3
EOF
    
    t5_1="; Ligand position restraints"
    t5_2="#ifdef POSRES"
    t5_3="#include \"posre_"$l".itp\""
    t5_4="#endif"
    
    sed -i "/; Include water topology/i$t5_1\n$t5_2\n$t5_3\n$t5_4\\"$'\n' topol.top
    echo " "
    echo "Gerando no indice o grupo Protein-"$L" e "$p_ion"_"$n_ion"_Water..."
    echo " "
    gmx make_ndx -f em.gro -o index.ndx <<EOF
    1 | 13
    14 | 15 | 17
    q
EOF

    #REALIZANDO O EQUILÍBRIO NVT:
    #CONFIRA AS CONFIGURAÇÕES DO ARQUIVO nvt.mdp
    echo "Executando a dinâmica de Equilíbrio a Temperatura constante NVT..."
    echo " "
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
    wait
    gmx mdrun -deffnm nvt -ntmpi 1 -ntomp "$n_threads" -nb gpu -bonded gpu -pme gpu -pmefft gpu -update gpu

    echo "Gerando o gráfico de Temperatura..."
    echo " "
    gmx energy -f nvt.edr -o temperature.xvg <<EOF
    16
    0
EOF

    #REALIZANDO O EQUILÍBRIO NPT:
    #CONFIRA AS CONFIGURAÇÕES DO ARQUIVO npt.mdp
    echo "Executando a dinâmica de Equilíbrio a Pressão constante NPT..."
    echo " "
    gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
    wait
    gmx mdrun -deffnm npt -ntmpi 1 -ntomp "$n_threads" -nb gpu -bonded gpu -pme gpu -pmefft gpu -update gpu

    echo "Gerando o gráfico de Pressão..."
    echo " "
    gmx energy -f npt.edr -o pressure.xvg <<EOF
    17
    0
EOF

    echo "Gerando o gráfico de Densidade..."
    echo " "
    gmx energy -f npt.edr -o density.xvg <<EOF
    23
    0
EOF

    #REALIZANDO A DINÂMICA - ETAPA DE PRODUÇÃO:
    echo "Executando o cálculo de Dinâmica Molecular..."
    echo " "
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_"$t_ns".tpr
    wait
    gmx mdrun -deffnm md_0_"$t_ns" -ntmpi 1 -ntomp "$n_threads" -nb gpu -bonded gpu -pme gpu -pmefft gpu -update gpu

    echo "Ações pós-dinâmica:"
    echo "Centralizando a proteína..."
    gmx trjconv -s md_0_"$t_ns".tpr -f md_0_"$t_ns".xtc -n index.ndx -o md_0_"$t_ns"_noPBC.xtc -pbc mol -center<<EOF
    1
    0
EOF

    echo "Calculando RMSD do complexo Poteína-Ligante..."
    gmx rms -s md_0_"$t_ns".tpr -f md_0_"$t_ns"_noPBC.xtc -n index.ndx -o rmsd_PL.xvg -tu ns<<EOF
    4
    19
EOF

    echo "Calculando DISTRMSD do complexo Poteína-Ligante..."
    gmx rmsdist -s md_0_"$t_ns".tpr -f md_0_"$t_ns"_noPBC.xtc -n index.ndx -o distrmsd.xvg<<EOF
    19
EOF

    echo "Calculando a flutuação dos resíduos da proteína..."
    gmx rmsf -s md_0_"$t_ns".tpr -f md_0_"$t_ns"_noPBC.xtc -n index.ndx -o rmsf.xvg -res<<EOF
    1
EOF

    echo "Calculando as ligações de hidrogênio entre proteína-ligante..."
    gmx hbond -f md_0_"$t_ns"_noPBC.xtc -s md_0_"$t_ns".tpr -n index.ndx -num hydrogen2-bond-intra-protein.xvg<<EOF
    1
    13
EOF

    echo "Calculando a área de superfície dos resíduos..."
    gmx sasa -f md_0_"$t_ns"_noPBC.xtc -s md_0_"$t_ns".tpr -o sas2.xvg -oa atomic2-sas.xvg -or residue2-sas.xvg<<EOF
    1
EOF

    echo "Calculando o centro de massa do complexo"
    gmx traj -f md_0_"$t_ns"_noPBC.xtc -s md_0_"$t_ns".tpr -n index.ndx -com yes -ox com_"$L"_PROT.xvg -x yes -y yes -z yes -nojump yes -tu ns<<EOF
    19
EOF
done
}
show_main_menu

##############################################################
# FUNÇÃO PARA EXECUTAR O CÁLCULO DE ENERGIA DE LIGAÇÃO (MMGBSA):
##############################################################
executa_mmgbsa() {
for lig in "$ligands"/*/; do
    L=$(basename "$lig")
    l="${L,,}"
    echo "Calculando MMGBSA..."
    AMBERHOME="home/moises/miniconda3"
    source $AMBERHOME/amber.sh
    conda activate gmxMMPBSA
    mpirun -np "$n_threads" gmx_MMPBSA MPI -O -i mmpbsa.in -cs md_0_"$t_ns".tpr -ci index.ndx -cg 1 13 -ct md_0_"$t_ns"_noPBC.xtc -cp topol.top &
    while true;
    do
        if pgrep gmx_MMPBSA_ana > /dev/null
        then
            pkill gmx_MMPBSA_ana
            break 
        else
            echo "O gmx_MMPBSA ainda não finalizou!"
            read -n 1 -t 600 -p "Pressione qualquer tecla para continuar..."
            date
        fi        
     done
    #É PRECISO EXECUTAR pkill gmx_MMPBSA_ana POIS, CASO CONTRÁRIO, ao final do cálculo de GBSA ele abrirá o MMPBSA_ana e ficará aguardando uma confirmação. Com isso o próximo ligante não iniciará pois o loop for para.
    echo "Término do cálculo de MMPBSA_MMGBSA..."
    # Fim do contador de tempo
    end_time=$(date +%s)

    # Cálculo do tempo decorrido
    end_time=$(date +%s)
    st=$((end_time - start_time))
    mt=$((st / 60))
    ht=$((mt / 60))
    s=$((st - mt * 60))
    m=$((mt - ht * 60))

    echo "Tempo decorrido: $st segundos" >> "$t"
    echo "Tempo decorrido: $ht horas:$m minutos:$s segundos" >> "$t"
done
}
show_main_menu
