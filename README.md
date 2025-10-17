# CrySPY-1.4.3_utility
- This is unofficial, but it has more examples than the original. I added them after confirming that it works.
- LAQA: only VASP, QE, or soiap for now

## Ternary system
- The "end_point" is a parameter that specifies the search range of the chemical potential. The standard is to obtain the DFT energy of a single molecule from databases such as Material Projects, literature, or your own calculations, and expand it by about +/- 1 eV.
- Fe-Pd-In
- Li-La-Ti
- Ni-Ga-Nd
- Fe-Pd-Pb
- Fe-Cr-Al
- Mn-Al-C
- Fe-Ni-Si
- Ti-Fe-Al

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
