make_input.py
-------------
This script generates `cryspy.in` and `POSCAR` files for evolutionary algorithm runs.
It supports 2 to 10 elements (or more) dynamically.

Usage:
    python3 make_input_lammps.py Element1 Element2 [Element3 ... Element10]

Example for a 5-element system:
    python3 make_input_lammps.py Al Ni Fe Co Zr


Note
-------------
- Do **NOT** delete the following files:
    Xx_tmp, calc_in, element_data_chgnet.txt, make_input.py

- Current implementation uses **1NN-MEAM** with:
    pair_coeff * * ./../../library_1nn.meam ${elem} NULL ${elem}
  This means only single-element parameters are applied.

- **Important for alloy systems:**
    - Add proper **MEAM parameter files** for cross-interactions.
    - Prepare alloy-specific entries in `library.meam`.
    - Consider switching to **2NN-MEAM** for better accuracy.

- **Critical limitation:**
    Any element **not defined in `library_1nn.meam` cannot be calculated**.
    Ensure all elements you specify have corresponding MEAM parameters in the library file.

- **Energy settings in `element_data_lammps.txt`:**
    The energy values are computed as:
        -esub + 0.043361254529175 * 2.5 * 2
    where `-esub` is taken from `library_1nn.meam`.
    This adjustment was chosen because:
        - It makes phase diagram plots easier to interpret visually.
        - The added term corresponds to ~25 kcal/mol, representing a thermal window similar to an oil bath heating condition where reactions typically occur.

- **Why this matters:**
    1NN-MEAM can provide a rough, qualitative phase diagram, which may be sufficient for initial screening and MD predictions.
    However, for reliable alloy behavior and phase stability, cross-element interactions and second-nearest neighbor terms are essential.

- **Future improvements:**
    - Validate results against DFT or experimental data.
    - Update MEAM libraries for multi-element systems.
    - Use 2NN-MEAM for more realistic thermodynamics and kinetics.