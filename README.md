# CrySPY-1.4.3_utility
- This is unofficial, but it has more examples than the original. I added them after confirming that it works.

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
