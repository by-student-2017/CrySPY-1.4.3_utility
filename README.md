# CrySPY-1.4.3_utility
- This is unofficial, but it has more examples than the original. I added them after confirming that it works.
- This is a calculation example. For practical use, it is recommended to increase the number of nat and tot_struc.
- LAQA: only VASP, QE, or soiap for now

## Note: Results and random numbers
- The random seed is changed so the results will not match those in the data. Therefore, even if the environment is almost the same and the user is the same, the results will be different.

## Calculated phase diagram
- The "end_point" is a parameter that specifies the search range of the chemical potential. The standard is to obtain the DFT energy of a single molecule from databases such as Material Projects, literature, or your own calculations, and expand it by about +/- 1 eV.

## How to use "plot3d.py", "plot3dn.py" or "pplot_v1.py"
- It is designed to correspond to a four-element system. The figure can be obtained by the following command.
```
python3 plot3d.py
```
```
python3 plot3dn.py
```
- The number of atoms of one element in the quaternary system is fixed (nat is fixed to a certain value) and the calculation is performed as a pseudo-ternary system, and the results are shown in the figures. The name "v1" is used because the guide lines have not yet been drawn properly. ("PTD" is an abbreviation for pseudo-ternary diagram.)
```
python3 pplot_v1.py
```

## Automatic generation of input files (EA-vc + CHGNet version)
- Since creating an input file manually can be cumbersome, a Python script called make_input.py is provided for CHGNet (see ase_chgnet_auto_EA-vc).
```
python3 make_input.py Al C Ni Fe
```
- Specify the elements after make_input.py. You can also create ternary or quinary systems.
- The file element_data_chgnet.txt is critical because it is used by ASE to generate POSCAR files.
- This file contains reference structures and lattice constants, which may include simplified or metastable phases for convenience (e.g., Boron as a diamond structure).
- These values are intended for efficient initial structure generation and may differ from actual ground-state phases.
- For accurate results, it is strongly recommended to re-optimize the final structures and energies using DFT or other first-principles methods.
- If necessary, you can edit element_data_chgnet.txt to replace entries with known stable phases or literature values.
- If the generated results appear incorrect, check and adjust this file accordingly.
- If you want to use B_aB12, just change B_dia to a directory named B_aB12, change POSCAR to correspond to the structure of aB12, and change the end_point in cryspy to the eV/atom of aB12. These can be done after "make_input.py". 

## Automatic generation of input files (EA-vc + DFTB+(GFN1-xTB) version)
- The "element_data_dftb.txt" has not been rewritten with the GFN1-xTB results in DFTB+ (I plan to fix this in the future if I have time). Therefore, an element_structure file will be created, so please calculate it and correct the end_point.

## How to use "repeat_cryspy.sh"
- "repeat_cryspy.sh" can be used anywhere. You can copy it and use it instead of cryspy as follows.
```
chmod +x repeat_cryspy.sh
bash repeat_cryspy.sh
```

## DFTB+ Installation
- Older versions of DFTB+ allow you to use Pressure with LBFGS, so in some cases you may want to consider using an older version. 
```
### DFTB+ ver. 24.1
cd $HOME
wget https://github.com/dftbplus/dftbplus/releases/download/24.1/dftbplus-24.1.x86_64-linux.tar.xz
tar -xvf dftbplus-24.1.x86_64-linux.tar.xz
echo 'export PATH=$PATH:$HOME/dftbplus-24.1.x86_64-linux/bin' >> ~/.bashrc
source ~/.bashrc
```

## Installation (CrySPY)
```
# Installation (CrySPY)
cd $HOME
pip install csp-cryspy==1.4.3
pip install numexpr==2.14.1
pip install bottleneck==1.6.0
pip install chgnet==0.4.2
mkdir work
cd work
git clone https://github.com/Tomoki-YAMASHITA/CrySPY_utility.git
```
## Note: Other method 1 (Tutorial files)
```
wget https://github.com/Tomoki-YAMASHITA/CrySPY_utility/archive/refs/heads/master.zip
sudo apt install unzip
unzip master.zip
mv CrySPY_utility-master CrySPY_utility
```
## Note: Other method 2 (Tutorial files) (This gihtub has more examples.)
```
wget https://github.com/by-student-2017/CrySPY-1.4.3_utility/archive/refs/heads/master.zip
sudo apt install unzip
unzip master.zip
mv CrySPY-1.4.3_utility-master CrySPY_utility
```
## For show svg fomat data
```
sudo apt update
sudo apt -y install eog
```

## Environments
```
[2025-10-16 16:02:16,307][cryspy_init][INFO]

Start CrySPY 1.4.3

[2025-10-16 16:02:16,307][cryspy_init][INFO] # ---------- Library version info
[2025-10-16 16:02:16,307][cryspy_init][INFO] ase version: 3.26.0
[2025-10-16 16:02:16,310][cryspy_init][INFO] pandas version: 2.3.3
[2025-10-16 16:02:16,311][cryspy_init][INFO] pymatgen version: 2025.10.7
[2025-10-16 16:02:16,311][cryspy_init][INFO] pyxtal version: 1.1.1
[2025-10-16 16:02:16,311][cryspy_init][INFO] # ---------- Read input file, cryspy.in
```
- BO (numba==0.62.1 case)
```
pip3 install physbo==2.2.0
pip3 install dscribe==2.1.2
pip install --upgrade numba
```
## Recalculation methods for low-cost calculations like Lammps (or CHGNet) and high-cost calculations like QE (VASP, DFTB+(xTB)).
- After performing the calculation using the low-cost method, copy "data/pkl_data/init_struct_data.pkl" and "data/pkl_data/opt_struc_data.pkl".
- Use the load_struc_flag = True option.
- The [option] load_struc_flag = True setting is only loaded the first time the program is run, so you don't need to worry about it from the second time onwards.

## The "kppvol" and actual k-points
- The "kppvol" is k-points/volume [A^-3].
- The actual k-points to be set can be determined using "kpt_check.py".
- A kppovl value of 40 is appropriate for rough calculations, while a value of 80 is appropriate for full calculations, so no special changes are necessary.

## References
- [1] CCMSハンズオン: CrySPY講習会 (2025/10/15)
- [2] ASEとCrySPYでサクッと結晶構造探索 #ASE - Qiita
- [3] Hitoshi GOMI's webpage
- [4] マテリアルズインフォマティクス 伊藤 聡(編) - 共立出版
- [5] 物性研スパコンohtakaで結晶構造探索を行う #Python - Qiita
- [6] CrySPY: a crystal structure prediction tool accelerated by machine learning
- [7] csp-cryspy · PyPI

---

## Phase diagrams attracting attention
### Ternary system
- Fe-Pd-In
- Li-La-Ti
- Ni-Ga-Nd
- Fe-Pd-Pb
- Fe-Cr-Al
- Mn-Al-C
- Fe-Ni-Si
- Ti-Fe-Al
---
- Li-Cu-Sb
- Zn-Sn-O
- Cu-Al-Mn
---
- Ti-V-Cr
- Ti-V-Al
- Ti-V-Ni
- Ti-V-H
- Mg-Ti-H
- Fe-Ti-H
---
- Ni-Cr-W
- Fe-Mo-B
- Zn-Sn-P
- Li-Zr-O
- Mg-Al-C
- Ti-Al-N
- Fe-Ni-Cr
---
- Be-Al-Mg
- Be-Ti-Zr
- Be-Ni-Cr
- Be-Li-Mg
- Be-Li-Al
- Be-Li-Ti
---
- Co-Cr-Al
- Ni-Co-Ti
- Fe-Co-Ni
- Mg-Zn-Ca
- Li-Ni-Mn
- Ti-Zr-Hf
- Cu-Zn-Al
- Nb-Ti-Al
- Si-Ge-Sn
- Zn-Ga-O
---
- Co-Ni-Al
- Ti-Nb-Zr
- Mo-Si-B
- Al-Ti-C
- Cr-Mn-Ni
- Ga-In-Zn
- Li-Si-P
- Mg-Sn-Ca
---
- La-Si-Mg
- Zr-Y-Ti
- Li-Zr-P
- Fe-Mo-B
---
- Li-Mg-Si
- Na-Li-Zr
- Ti-Nb-Mo
- Cr-Al-Si
- Ga-Zn-Sn
- Nb-Mo-Si
- Sc-Y-Zr
---
- Ca-Mg-Al
- Ti-Si-C
- Al-Cr-Nb
- Li-Mg-Zn
- Fe-Al-Si
- Ga-Zn-Sn
---
- Al-Nb-Mo
- Ti-Cr-Mo
- Ca-Mg-Zn
---
- Nb-Mo-Zr
- Nb-Zr-Hf
- Nb-Zr-Al
- Nb-Zr-Si
---
- Sc-Y-Zr
- Sc-Y-Ti
- Sc-Y-Al
- Sc-Y-Nb
---
- Ti-Zr-Hf
- Ti-Nb-Hf
- Hf-Zr-Si

### Quaternary system
- Fe-Cr-Al-Y
- Mn-Al-C-Ni
- Li-La-Ti-Zr
- Fe-Pd-In-Sn
---
- Li-La-Zr-O
- Fe-Ni-Cr-Al
- Ti-Al-Nb-Mo
---
- Co-Cr-Fe-Ni
- Al-Co-Cr-Ni
- Li-Ni-Mn-Co
- Ti-Zr-Hf-Nb
- Fe-Co-Ni-Cr
- Cu-Al-Mn-Ni
- Zn-Ga-In-O
- Si-Ge-Sn-Pb
- Mg-Zn-Ca-Sr
---
- Co-Cr-Ni-Mo
- Al-Ti-Nb-Zr
- Li-Na-K-Cl
- Ga-In-Zn-O
- Fe-Ni-Co-Cr
- Mo-Si-B-Ti
---
- Li-Mn-Co-Ni
- Ga-Zn-In-O
----
- Li-Mg-Si-O
- Ti-Nb-Mo-Zr
- Co-Cr-Ni-Mo
- Ga-Zn-Sn-O
- Li-Na-K-F
- Mo-Nb-Si-B
---
- Al-Cr-Nb-Ti
- Li-Mg-Si-O
- Fe-Al-Cr-Si
- Ti-Nb-Mo-Zr
---
- Ga-Au-Pt-Gd
- Li-La-Zr-O
- Zn-Ga-In-O
- Mo-Nb-Si-B
- Li-Na-K-Cl
---
- Li-Mg-Al-O
- Mo-Nb-Ti-Zr
- Sc-Y-Ti-Zr
- Hf-Nb-Mo-Zr
