make_input.py
-------------
This script generates cryspy.in and POSCAR files for evolutionary algorithm runs.
It supports 2 to 10 elements (or more) dynamically.

Usage:
    python3 make_input_qe.py Element1 Element2 [Element3 ... Element10]

Example for a 2-element system:
    python3 make_input_qe.py Fe C


Note
-------------
Do not delete "Xx_tmp, calc_in, element_data_chgnet.txt, make_input_chgnet.py".
The "eV/atom" is the enthalpy. (QE does not yet support reading enthalpy using the energy_step_flag option.)
Enthalpy: "pv_term = True" and "press = 400". see "qe_Sr4O4_RS_pv_term"
The unit of pressure in QE is [kbar] (1 [bar] = 0.1 [MPa]).
