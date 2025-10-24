make_input.py
-------------
This script generates cryspy.in and POSCAR files for evolutionary algorithm runs.
It supports 2 to 10 elements (or more) dynamically.

Usage:
    bash run.sh Element1 Element2 [Element3 ... Element10]

Example for a 2-element system:
    bash run.sh Fe C


Note
-------------
Do not delete "Xx_tmp, calc_in, element_data_chgnet.txt, make_input_chgnet.py".

CHGNet and MACE are more compatible with VASP than with Quantum Espresso, 
so it would be better to rewrite the code for VASP. However, I'm just a low-level user, so I can't use VASP.