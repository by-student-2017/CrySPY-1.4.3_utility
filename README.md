# CrySPY-1.4.3_utility
- This is unofficial, but it has more examples than the original. I added them after confirming that it works.
- LAQA: only VASP, QE, or soiap for now

## Calculated phase diagram
- The "end_point" is a parameter that specifies the search range of the chemical potential. The standard is to obtain the DFT energy of a single molecule from databases such as Material Projects, literature, or your own calculations, and expand it by about +/- 1 eV.

## How to use "plot3d.py" or "plot3dn.py"
- It is designed to correspond to a four-element system. The figure can be obtained by the following command.
```
python3 plot3d.py
```
```
python3 plot3dn.py
```

# How to use "repeat_cryspy.sh"
- "repeat_cryspy.sh" can be used anywhere. You can copy it and use it instead of cryspy as follows.
```
chmod +x repeat_cryspy.sh
bash repeat_cryspy.sh
```

# DFTB+ Installation
- Older versions of DFTB+ allow you to use Pressure with LBFGS, so in some cases you may want to consider using an older version. 
```
### DFTB+ ver. 24.1
cd $HOME
wget https://github.com/dftbplus/dftbplus/releases/download/24.1/dftbplus-24.1.x86_64-linux.tar.xz
tar -xvf dftbplus-24.1.x86_64-linux.tar.xz
echo 'export PATH=$PATH:$HOME/dftbplus-24.1.x86_64-linux/bin' >> ~/.bashrc
source ~/.bashrc
```

# Environments
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

## Ternary system
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

## Quaternary system
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
