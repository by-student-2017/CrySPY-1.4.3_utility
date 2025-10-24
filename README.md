# CrySPY-1.4.3_utility
- This is unofficial, but it has more examples than the original. I added them after confirming that it works.
- This is a calculation example. For practical use, it is recommended to increase the number of nat and tot_struc.
- LAQA: only VASP, QE, or soiap for now

## Results and random numbers
- The random seed is changed so the results will not match those in the data. Therefore, even if the environment is almost the same and the user is the same, the results will be different.

## Calculated phase diagram
- The "end_point" is a parameter that specifies the search range of the chemical potential. The standard is to obtain the DFT energy of a single molecule from databases such as Material Projects, literature, or your own calculations, and expand it by about +/- 1 eV.
- The element_structure file in the example is not actually necessary for EA-vc. It is output for reference. Therefore, the make_input_*.py shown below can be made simpler. Since simpler is often better, users may wish to simplify it themselves.
- The general rule for state diagrams is to arrange them in alphabetical order. General names may not always be arranged in alphabetical order. Just be aware of the general rule.
- Temperature control is performed using molecular dynamics (MD). The effects of lattice vibrations (phonons) are taken into account, but the calculation takes time.
- For pressure application at 0 K, it is also a good idea to use QE or the old DFTB+. It is also a good method to find candidates from there and examine them in more detail using MD etc.
- When applying temperature and pressure, the diagram will only be accurate if the energy and end_point values in element_data_*.txt are adjusted to reflect the crystal structure at those conditions. If your goal is simply to highlight differences in stable structures, you can make the diagram easier to read by setting these values to zero or to a slightly positive number.

## How to use "plot3d.py", "plot3dn.py", "pplot_v1.py", or "plot1d.py"
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
- Use "cryspy_rslt" obtained in any calculation to plot a diagram of Structure ID vs. Energy (eV/atom).
```
python3 plot1d.py
```
- Instead of plot1d.py, a gnuplot version is also available. The usage is as follows (-persist can be omitted):
```
gnuplot -persist plot1d.gpl
```
- Gnuplot can also be used to create ternary and quaternary phase diagrams, but changing the axis angles is extremely tedious, so it is not provided here. (Readers with ample time might want to try making them, although researchers with such time are probably a rare breed these days.)

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
- GFN2-xTB can also be specified in dftb_in.hsd, but GFN2-xTB is difficult to handle due to its computational cost and convergence. In the paper, good results were obtained with GFN1-xTB for MOF, so here we create an example using GFN1-xTB.
- Interested readers may want to try using GFN2-xTB or the skf file, etc.

## Automatic generation of input files (EA-vc + QE version)
- The "element_data_qe.txt" does not reflect the results of QE (we plan to fix this in the future if we have time). Therefore, the element_structure file will be created, so please perform the calculation and correct the end_point (element_structure may also not be created correctly, so please check).
- QE and VASP can also be discussed in terms of enthalpy by entering only pressure. Please refer to the official manual.

## Automatic generation of input files (EA-vc + Lammps version)
- The Lammps results are not reflected in "element_data_lammps.txt" (we plan to fix this in the future if we have time). Therefore, the element_structure file will be created, so please run the calculation and correct the end_point (the element_structure may not have been created correctly, so please check).
- The Lammps EAM potential code has been rewritten to fully match the Zhou 2004 data published by NIST, and parameters for previously reported elements have been added. It is theoretically more expressive than Lennard-Jones (LJ), so it may be worth considering using it. However, not all elements are included.

## Automatic generation of input files (EA-vc + MACE version)
- The MACE results are not reflected in "element_data_mace.txt" (I plan to fix this in the future if I have time). Therefore, the element_structure file will be created, so please run the calculation and correct the end_point (it's possible that the element_structure was not created correctly, so please check). As with CHGNet, I think it's fine to match it to the VASP results. It's also a good idea to organize it with Material Projects data.

## Enthalpy (pv_term) (such as in high-pressure research)
- Examples of input files are located in the directory labeled "pv_term".
- Currently (24/Oct/2025), only QE and VASP are supported. QE does not yet support reading enthalpy using the energy_step_flag option. (1 [kbar] = 0.1 [GPa])

## Random Search (RS)
- Even after a series of calculations has finished, you can continue the calculations by increasing the value of tot_struc in cryspy.in.
- This is an important method if you want to find the most stable phase rather than a phase diagram.

## How to use "repeat_cryspy.sh"
- "repeat_cryspy.sh" can be used anywhere. You can copy it and use it instead of cryspy as follows.
```
chmod +x repeat_cryspy.sh
bash repeat_cryspy.sh
```

## Linking different calculation codes
- We built a system that allows calculations to continue using different calculation codes by utilizing "load_struc_flag=True".
- I'm using Quantum ESPRESSO (QE) because I don't have a VASP license or the budget for a high-performance computer.
Since CHGNet and MACE utilize VASP data, calculations don't proceed as smoothly with QE (whereas structural optimization would likely be much faster with VASP).
- I think that by improving "ase_chgnet_qe_auto_EA-vc", it would be relatively easy to create combinations with other calculation codes, so I hope users will do their best with VASP on their own.

## Useful features
- Deleting old files: "cryspy -c" or "cryspy --clean"
- Enthalpy: "pv_term = True" and "press = 400". see "qe_Sr4O4_RS_pv_term"
- Replot: "cryspy-Eplot" or "cryspy-calc-convex-hull 1"
- Jupyter Notebook: https://tomoki-yamashita.github.io/CrySPY_doc/tutorial/interactive/index.html

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
### Usage eog
```
eog data/convex_hull/conv_hull_gen_1.svg
```
## For gnuplot (plot1d.gpl)
```
sudo apt -y install gnuplot
```
## MACE Installation
```
pip install mace-torch==0.3.14
pip install mace-models==0.1.6
```
### MACE (GPU case)
```
pip install cuequivariance cuequivariance-torch
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
- BO (numba==0.62.1 case) (using COMBO code)
```
pip3 install physbo==2.2.0
pip3 install dscribe==2.1.2
pip install --upgrade numba
```
### PC specs used for test
- OS: Microsoft Windows 11 Home 64 bit
- BIOS: 1.14.0
- CPU： 12th Gen Intel(R) Core(TM) i7-12700
- Base Board：0R6PCT (A01)
- Memory：32 GB
- GPU: NVIDIA GeForce RTX3070
- WSL2: VERSION="22.04.1 LTS (Jammy Jellyfish)"
- Python 3.10.12
### requirements.txt
```
pip freeze > requirements.txt
```

## The "kppvol" and actual k-points
- The "kppvol" is k-points/volume [A^-3].
- The actual k-points to be set can be determined using "kpt_check.py".
- A kppovl value of 40 is appropriate for rough calculations, while a value of 80 is appropriate for full calculations, so no special changes are necessary.

## kppvol (k-points per volume) Standards and Accuracy Guidelines
- Energy convergence depends on the system size and band structure, but as a rule of thumb:
- kppvol = 30–50: Accuracy of a few meV/atom in solids
- kppvol = 60–80: Often aiming for 1 meV/atom or less
- (kppvol = 25: 5-10 meV/atom, kppvol = 35: 2-5 meV/atom, kppvol = 40: 2-4 meV/atom)
### In metallic systems, introducing smearing (Methfessel-Paxton, Gaussian, etc.) reduces the energy oscillations caused by discontinuities in the Fermi surface, making it easier to maintain accuracy even with a somewhat coarse k-point mesh.
- Metals with smearing: kppvol 35–45, approximately a few meV/atom
- High accuracy (phase transitions and small energy differences): kppvol ≥ 60 recommended
- kppvol = 35 is somewhat coarse, and typically results in an error of 2–5 meV/atom.
### Changing the calculation conditions for the purpose of exploration
- For magnetic systems, the energy difference is very small (often less than a few meV/atom), so the accuracy requirements are high. Especially when comparing magnetic order (FM vs AFM) or the stability of spin structures, an error of 1 meV/atom level can have a significant impact on the results. However, the accuracy requirements vary depending on your application in CrySPY:
- Initial structure search (EA or random search): A rough search is OK. kppvol = 30–40 is sufficient. Since the goal is candidate generation, some error is acceptable.
- Final evaluation (determination of magnetic order and relative stability): High accuracy is required. kppvol >= 60 is recommended. Magnetic systems have complex band structures, and smearing has limited effectiveness, so a dense k-point is required.

## References
- [1] CCMSハンズオン: CrySPY講習会 (2025/10/15)
- [2] マテリアルズインフォマティクス 伊藤 聡(編) - 共立出版
- [3] Hitoshi GOMI's webpage
- [4] 物性研スパコンohtakaで結晶構造探索を行う #Python - Qiita
- [5] ASEとCrySPYでサクッと結晶構造探索 #ASE - Qiita
- [6] CrySPY: a crystal structure prediction tool accelerated by machine learning
- [7] csp-cryspy · PyPI

---

## Phase Diagrams Attracting Attention

Phase diagrams are essential for understanding phase stability and material properties. The following points highlight current challenges and recommendations:

- **International Collaboration:**  
  Phase diagrams for precious metals, rare earths, and other elements that are difficult to obtain experimentally should ideally be compiled at multiple computational levels (e.g., CALPHAD, first-principles) through cooperation between countries.

- **Organization Rule:**  
  The general convention for arranging phase diagrams is alphabetical order. While common names may occasionally deviate, awareness of this rule ensures consistency.

- **Current Status and Future Needs:**  
  Phase diagrams are being collected by organizations such as **NIMS**, but the available information is incomplete. It would be highly beneficial to compile these datasets using integrated platforms such as **Phase0** or similar frameworks to enable collaborative development and global accessibility.

- **Next Research Stage:**  
  Once phase diagrams are clarified, researchers can advance to the next stage: performing molecular dynamics (MD) simulations in composition regions of interest. These simulations enable the study of dislocation behavior, mechanical properties, and environmental effects, providing deeper insights into material performance under realistic conditions.  
  **Note:** CrySPY focuses on crystal structure prediction. For MD and related property evaluations, tools such as **ASE**, **LAMMPS**, **Quantum ESPRESSO (QE)**, **VASP**, **Siesta**, **CP2K**, and **OpenMX** should be utilized.
  + #### Related Topics
  + Lammps: https://github.com/by-student-2017/lammps_education_metal_win
  + Lammps (+ charge): https://github.com/by-student-2017/lammps_education_reaxff_win
  + DFTB+: https://github.com/by-student-2017/DFTBplus-v.23.1-examples
    + [DFTB.org Parameter Repository:「Periodic Table Baseline Parameter Set (PTBP)」](https://dftb.org/parameters/download.html): CrySPY -> CHGNet, MACE or DFTB+ (fast search with GFN1-xTB) -> Evaluate top candidates with DFTB+ (PTBP) for electric field properties -> Finally evaluate with DFT.
    + SOAC: [GitHub - Voss-Lab/SK_repository](https://github.com/Voss-Lab/SK_repository/tree/main)
    + Handbook of the Band Structure of Elemental Solids
    + SlaKoNet: A Unified Slater-Koster Tight-Binding Framework: https://github.com/atomgptlab/slakonet
  + CP2k: https://github.com/by-student-2017/CP2k-v.9.1-examples
  + CALPHAD: https://github.com/by-student-2017/OpenCALPHAD-v6-examples
  + CALPHAD: https://github.com/by-student-2017/pycalphad-v.0.10.3-examples
  + Phase field: https://github.com/by-student-2017/Phase-field

# Ternary Alloy Systems

This document lists selected ternary alloy systems of interest for advanced materials research. Each entry includes a short description and the current status of phase diagram development.

---

### Fe-Pd-In
Intermetallic system with potential catalytic and electronic properties.  
*Phase Diagram:* Not fully established.

### Li-La-Ti
Lithium-based system for solid-state electrolytes and advanced ceramics.  
*Phase Diagram:* Limited data.

### Ni-Ga-Nd
System with possible magnetic and electronic applications.  
*Phase Diagram:* Sparse.

### Fe-Pd-Pb
Intermetallic system with potential for corrosion-resistant coatings.  
*Phase Diagram:* Unknown.

### Fe-Cr-Al
Oxidation-resistant alloy system for high-temperature applications.  
*Phase Diagram:* Well-studied in binary subsets; ternary incomplete.

### Mn-Al-C
Lightweight alloy system with possible magnetic properties.  
*Phase Diagram:* Developing.

### Fe-Ni-Si
System relevant to steels and magnetic materials.  
*Phase Diagram:* Partial data.

### Ti-Fe-Al
Refractory alloy system for structural applications.  
*Phase Diagram:* Sparse.

---

### Li-Cu-Sb
Lithium alloy system for battery anodes.  
*Phase Diagram:* Limited.

### Zn-Sn-O
Transparent oxide semiconductor system.  
*Phase Diagram:* Developing.

### Cu-Al-Mn
Shape-memory alloy system for actuators.  
*Phase Diagram:* Partial.

---

### Ti-V-Cr
Refractory alloy system for high-temperature uses.  
*Phase Diagram:* Sparse.

### Ti-V-Al
Lightweight alloy system for aerospace.  
*Phase Diagram:* Limited.

### Ti-V-Ni
System with hydrogen storage potential.  
*Phase Diagram:* Developing.

### Ti-V-H
Hydride-forming system for energy storage.  
*Phase Diagram:* Sparse.

### Mg-Ti-H
Hydride system for hydrogen storage.  
*Phase Diagram:* Limited.

### Fe-Ti-H
Hydride-forming alloy system.  
*Phase Diagram:* Sparse.

---

### Ni-Cr-W
High-temperature alloy system for corrosion resistance.  
*Phase Diagram:* Partial.

### Fe-Mo-B
System for hard coatings and wear resistance.  
*Phase Diagram:* Sparse.

### Zn-Sn-P
Semiconductor alloy system for electronics.  
*Phase Diagram:* Developing.

### Li-Zr-O
Lithium zirconate system for solid electrolytes.  
*Phase Diagram:* Limited.

### Mg-Al-C
Lightweight alloy system for structural applications.  
*Phase Diagram:* Sparse.

### Ti-Al-N
System for hard coatings and cutting tools.  
*Phase Diagram:* Well-studied in binary; ternary incomplete.

### Fe-Ni-Cr
Basis for stainless steels and superalloys.  
*Phase Diagram:* Well-characterized.

---

### Be-Al-Mg
Lightweight alloy system for aerospace.  
*Phase Diagram:* Sparse.

### Be-Ti-Zr
Refractory alloy system for extreme environments.  
*Phase Diagram:* Limited.

### Be-Ni-Cr
System for corrosion-resistant alloys.  
*Phase Diagram:* Sparse.

### Be-Li-Mg
Lightweight alloy system for structural uses.  
*Phase Diagram:* Limited.

### Be-Li-Al
Lightweight alloy system for aerospace.  
*Phase Diagram:* Sparse.

### Be-Li-Ti
System for advanced lightweight alloys.  
*Phase Diagram:* Sparse.

---

### Co-Cr-Al
Oxidation-resistant alloy system for coatings.  
*Phase Diagram:* Developing.

### Ni-Co-Ti
System for high-temperature alloys.  
*Phase Diagram:* Sparse.

### Fe-Co-Ni
Magnetic alloy system for soft magnets.  
*Phase Diagram:* Well-studied.

### Mg-Zn-Ca
Biodegradable alloy system for medical implants.  
*Phase Diagram:* Limited.

### Li-Ni-Mn
Cathode material system for Li-ion batteries.  
*Phase Diagram:* Well-characterized.

### Ti-Zr-Hf
Refractory alloy system for high-temperature uses.  
*Phase Diagram:* Sparse.

### Cu-Zn-Al
Shape-memory alloy system for actuators.  
*Phase Diagram:* Partial.

### Nb-Ti-Al
Refractory alloy system for aerospace.  
*Phase Diagram:* Sparse.

### Si-Ge-Sn
Semiconductor alloy system for thermoelectrics.  
*Phase Diagram:* Developing.

### Zn-Ga-O
Transparent oxide semiconductor system.  
*Phase Diagram:* Limited.

---

### Co-Ni-Al
System for high-temperature alloys.  
*Phase Diagram:* Sparse.

### Ti-Nb-Zr
Refractory alloy system for biomedical and aerospace.  
*Phase Diagram:* Limited.

### Mo-Si-B
System for oxidation-resistant coatings.  
*Phase Diagram:* Sparse.

### Al-Ti-C
System for hard coatings and cutting tools.  
*Phase Diagram:* Developing.

### Cr-Mn-Ni
System relevant to stainless steels.  
*Phase Diagram:* Well-characterized.

### Ga-In-Zn
Transparent oxide semiconductor system.  
*Phase Diagram:* Developing.

### Li-Si-P
System for battery anodes and solid electrolytes.  
*Phase Diagram:* Limited.

### Mg-Sn-Ca
Biodegradable alloy system for implants.  
*Phase Diagram:* Sparse.

---

### La-Si-Mg
System for advanced ceramics.  
*Phase Diagram:* Limited.

### Zr-Y-Ti
System for high-temperature ceramics and alloys.  
*Phase Diagram:* Sparse.

### Li-Zr-P
System for solid electrolytes.  
*Phase Diagram:* Limited.

### Fe-Mo-B
System for hard coatings and wear resistance.  
*Phase Diagram:* Sparse.

---

### Li-Mg-Si
System for lightweight alloys and battery materials.  
*Phase Diagram:* Limited.

### Na-Li-Zr
System for solid electrolytes.  
*Phase Diagram:* Sparse.

### Ti-Nb-Mo
Refractory alloy system for high-temperature uses.  
*Phase Diagram:* Sparse.

### Cr-Al-Si
System for oxidation-resistant coatings.  
*Phase Diagram:* Developing.

### Ga-Zn-Sn
Transparent oxide semiconductor system.  
*Phase Diagram:* Limited.

### Nb-Mo-Si
System for refractory alloys.  
*Phase Diagram:* Sparse.

### Sc-Y-Zr
Rare-earth alloy system for advanced applications.  
*Phase Diagram:* Limited.

---

### Ca-Mg-Al
Lightweight alloy system for structural uses.  
*Phase Diagram:* Sparse.

### Ti-Si-C
System for hard coatings and ceramics.  
*Phase Diagram:* Developing.

### Al-Cr-Nb
Refractory alloy system for aerospace.  
*Phase Diagram:* Sparse.

### Li-Mg-Zn
Lightweight alloy system for structural uses.  
*Phase Diagram:* Limited.

### Fe-Al-Si
System for oxidation-resistant alloys.  
*Phase Diagram:* Developing.

### Ga-Zn-Sn
Transparent oxide semiconductor system.  
*Phase Diagram:* Limited.

---

### Al-Nb-Mo
Refractory alloy system for high-temperature uses.  
*Phase Diagram:* Sparse.

### Ti-Cr-Mo
System for corrosion-resistant alloys.  
*Phase Diagram:* Limited.

### Ca-Mg-Zn
Biodegradable alloy system for implants.  
*Phase Diagram:* Sparse.

---

### Nb-Mo-Zr
Refractory alloy system for high-temperature uses.  
*Phase Diagram:* Sparse.

### Nb-Zr-Hf
Refractory alloy system for extreme environments.  
*Phase Diagram:* Limited.

### Nb-Zr-Al
System for advanced structural alloys.  
*Phase Diagram:* Sparse.

### Nb-Zr-Si
System for oxidation-resistant coatings.  
*Phase Diagram:* Limited.

---

### Sc-Y-Zr
Rare-earth alloy system for advanced applications.  
*Phase Diagram:* Limited.

### Sc-Y-Ti
System for advanced lightweight alloys.  
*Phase Diagram:* Sparse.

### Sc-Y-Al
System for advanced structural alloys.  
*Phase Diagram:* Limited.

### Sc-Y-Nb
System for refractory alloys.  
*Phase Diagram:* Sparse.

---

### Ti-Zr-Hf
Refractory alloy system for high-temperature uses.  
*Phase Diagram:* Sparse.

### Ti-Nb-Hf
System for biomedical and aerospace alloys.  
*Phase Diagram:* Limited.

### Hf-Zr-Si
System for oxidation-resistant coatings.  
*Phase Diagram:* Sparse.

# Quaternary Alloy Systems

This document lists selected quaternary alloy systems of interest for advanced materials research. Each entry includes a short description and the current status of phase diagram development.

---

### Fe-Cr-Al-Y
High-temperature alloy system with excellent oxidation resistance.  
*Phase Diagram:* Not fully established.

### Mn-Al-C-Ni
Candidate for lightweight structural materials; magnetic properties possible.  
*Phase Diagram:* Incomplete.

### Li-La-Ti-Zr
Lithium-based system for solid-state electrolytes and advanced ceramics.  
*Phase Diagram:* Limited data.

### Fe-Pd-In-Sn
Potential for intermetallic compounds with catalytic properties.  
*Phase Diagram:* Not well documented.

---

### Li-La-Zr-O
Known for lithium-conducting ceramics (e.g., LLZO).  
*Phase Diagram:* Partially studied.

### Fe-Ni-Cr-Al
Basis for superalloys used in high-temperature applications.  
*Phase Diagram:* Well-studied for ternary subsets; quaternary incomplete.

### Ti-Al-Nb-Mo
Refractory alloy system for aerospace and turbine blades.  
*Phase Diagram:* Sparse data.

---

### Co-Cr-Fe-Ni
Core of high-entropy alloys with excellent mechanical properties.  
*Phase Diagram:* Under development.

### Al-Co-Cr-Ni
High-entropy alloy system with oxidation resistance.  
*Phase Diagram:* Incomplete.

### Li-Ni-Mn-Co
Cathode material system for Li-ion batteries (NMC).  
*Phase Diagram:* Well-characterized for practical compositions.

### Ti-Zr-Hf-Nb
Refractory HEA system for extreme environments.  
*Phase Diagram:* Limited experimental data.

### Fe-Co-Ni-Cr
High-entropy alloy system with balanced strength and ductility.  
*Phase Diagram:* Developing.

### Cu-Al-Mn-Ni
Shape-memory alloy system for actuators.  
*Phase Diagram:* Partial data.

### Zn-Ga-In-O
Transparent conductive oxide system for electronics.  
*Phase Diagram:* Not fully mapped.

### Si-Ge-Sn-Pb
Semiconductor alloy system for thermoelectric applications.  
*Phase Diagram:* Sparse.

### Mg-Zn-Ca-Sr
Biodegradable alloy system for medical implants.  
*Phase Diagram:* Limited.

---

### Co-Cr-Ni-Mo
Corrosion-resistant alloy system for biomedical and industrial uses.  
*Phase Diagram:* Incomplete.

### Al-Ti-Nb-Zr
High-temperature alloy system for aerospace.  
*Phase Diagram:* Sparse.

### Li-Na-K-Cl
Alkali halide system for molten salt applications.  
*Phase Diagram:* Well-known for binary/ternary; quaternary less studied.

### Ga-In-Zn-O
Transparent oxide semiconductor system.  
*Phase Diagram:* Developing.

### Fe-Ni-Co-Cr
High-entropy alloy system for structural applications.  
*Phase Diagram:* Partial.

### Mo-Si-B-Ti
Refractory alloy system for oxidation-resistant coatings.  
*Phase Diagram:* Sparse.

---

### Li-Mn-Co-Ni
Cathode material system for Li-ion batteries.  
*Phase Diagram:* Well-characterized for practical ranges.

### Ga-Zn-In-O
Transparent conductive oxide system.  
*Phase Diagram:* Limited.

---

### Li-Mg-Si-O
Lithium silicate-based system for solid electrolytes.  
*Phase Diagram:* Incomplete.

### Ti-Nb-Mo-Zr
Refractory alloy system for high-temperature structural uses.  
*Phase Diagram:* Sparse.

### Ga-Zn-Sn-O
Transparent oxide semiconductor system.  
*Phase Diagram:* Developing.

### Li-Na-K-F
Alkali fluoride system for molten salts.  
*Phase Diagram:* Limited.

### Mo-Nb-Si-B
Refractory alloy system for oxidation-resistant coatings.  
*Phase Diagram:* Sparse.

---

### Al-Cr-Nb-Ti
High-temperature alloy system for aerospace.  
*Phase Diagram:* Sparse.

### Fe-Al-Cr-Si
Oxidation-resistant alloy system for structural uses.  
*Phase Diagram:* Limited.

---

### Ga-Au-Pt-Gd
Exotic alloy system with potential magnetic and catalytic properties.  
*Phase Diagram:* Unknown.

### Zn-Ga-In-O
Transparent conductive oxide system for electronics.  
*Phase Diagram:* Not fully mapped.

---

### Li-Mg-Al-O
Lithium-aluminate-based system for ceramics and electrolytes.  
*Phase Diagram:* Limited.

### Mo-Nb-Ti-Zr
Refractory alloy system for high-temperature structural uses.  
*Phase Diagram:* Sparse.

### Sc-Y-Ti-Zr
Rare-earth and transition metal system for advanced alloys.  
*Phase Diagram:* Incomplete.

### Hf-Nb-Mo-Zr
Refractory alloy system for extreme environments.  
*Phase Diagram:* Not yet fully established.

---

# Undeveloped
- SevenNet (It's difficult for me.)
```
pip install e3nn>=0.5.0
pip install sevenn==0.11.2
```
- ORB v3 (It's difficult for me.)
```
pip install orb-models==0.5.5
pip install --extra-index-url=https://pypi.nvidia.com "cuml-cu12==25.2.*"
```
- Supports crystal, molecular, and surface systems
```
pip install matgl
```
