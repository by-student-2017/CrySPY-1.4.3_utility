# CrySPY-1.4.3_utility
- This is unofficial, but it has more examples than the original. I added them after confirming that it works.

- Operating environment
```
[2025-10-16 16:02:16,307][cryspy_init][INFO]


Start CrySPY 1.4.3


[2025-10-16 16:02:16,307][cryspy_init][INFO] # ---------- Library version info
[2025-10-16 16:02:16,307][cryspy_init][INFO] ase version: 3.26.0
[2025-10-16 16:02:16,310][cryspy_init][INFO] pandas version: 2.3.3
[2025-10-16 16:02:16,311][cryspy_init][INFO] pymatgen version: 2025.10.7
[2025-10-16 16:02:16,311][cryspy_init][INFO] pyxtal version: 1.1.1
[2025-10-16 16:02:16,311][cryspy_init][INFO] # ---------- Read input file, cryspy.in
[2025-10-16 16:02:16,312][write_input][INFO] [basic]
[2025-10-16 16:02:16,312][write_input][INFO] algo = EA-vc
[2025-10-16 16:02:16,312][write_input][INFO] calc_code = LAMMPS
[2025-10-16 16:02:16,312][write_input][INFO] nstage = 1
[2025-10-16 16:02:16,313][write_input][INFO] njob = 5
[2025-10-16 16:02:16,313][write_input][INFO] jobcmd = bash
[2025-10-16 16:02:16,313][write_input][INFO] jobfile = job_cryspy
[2025-10-16 16:02:16,313][write_input][INFO]
[2025-10-16 16:02:16,313][write_input][INFO] [structure]
[2025-10-16 16:02:16,313][write_input][INFO] struc_mode = crystal
[2025-10-16 16:02:16,313][write_input][INFO] atype = ('Cu', 'Ag', 'Au')
[2025-10-16 16:02:16,313][write_input][INFO] mindist_factor = 1.0
[2025-10-16 16:02:16,313][write_input][INFO] vol_factor = 1.1
[2025-10-16 16:02:16,313][write_input][INFO] symprec = 0.01
[2025-10-16 16:02:16,313][write_input][INFO] spgnum = all
[2025-10-16 16:02:16,313][write_input][INFO] use_find_wy = False
[2025-10-16 16:02:16,313][write_input][INFO] ll_nat = (0, 0, 0)
[2025-10-16 16:02:16,313][write_input][INFO] ul_nat = (8, 8, 8)
[2025-10-16 16:02:16,313][write_input][INFO]
[2025-10-16 16:02:16,313][write_input][INFO] [option]
[2025-10-16 16:02:16,313][write_input][INFO] check_mindist_opt = True
[2025-10-16 16:02:16,313][write_input][INFO] stop_chkpt = 0
[2025-10-16 16:02:16,313][write_input][INFO] load_struc_flag = False
[2025-10-16 16:02:16,313][write_input][INFO] stop_next_struc = False
[2025-10-16 16:02:16,313][write_input][INFO] append_struc_ea = False
[2025-10-16 16:02:16,313][write_input][INFO] energy_step_flag = False
[2025-10-16 16:02:16,313][write_input][INFO] struc_step_flag = False
[2025-10-16 16:02:16,313][write_input][INFO] force_step_flag = False
[2025-10-16 16:02:16,313][write_input][INFO] stress_step_flag = False
[2025-10-16 16:02:16,313][write_input][INFO]
[2025-10-16 16:02:16,313][write_input][INFO] [EA]
[2025-10-16 16:02:16,314][write_input][INFO] n_pop = 20
[2025-10-16 16:02:16,314][write_input][INFO] n_crsov = 5
[2025-10-16 16:02:16,314][write_input][INFO] n_perm = 2
[2025-10-16 16:02:16,314][write_input][INFO] n_strain = 2
[2025-10-16 16:02:16,314][write_input][INFO] n_rand = 2
[2025-10-16 16:02:16,314][write_input][INFO] n_elite = 2
[2025-10-16 16:02:16,314][write_input][INFO] fit_reverse = False
[2025-10-16 16:02:16,314][write_input][INFO] n_fittest = 10
[2025-10-16 16:02:16,314][write_input][INFO] slct_func = TNM
[2025-10-16 16:02:16,314][write_input][INFO] t_size = 2
[2025-10-16 16:02:16,314][write_input][INFO] crs_lat = random
[2025-10-16 16:02:16,314][write_input][INFO] nat_diff_tole = 4
[2025-10-16 16:02:16,314][write_input][INFO] ntimes = 1
[2025-10-16 16:02:16,314][write_input][INFO] sigma_st = 0.5
[2025-10-16 16:02:16,314][write_input][INFO] maxcnt_ea = 50
[2025-10-16 16:02:16,314][write_input][INFO] maxgen_ea = 0
[2025-10-16 16:02:16,314][write_input][INFO] n_add = 3
[2025-10-16 16:02:16,314][write_input][INFO] add_max = 3
[2025-10-16 16:02:16,314][write_input][INFO] n_elim = 3
[2025-10-16 16:02:16,314][write_input][INFO] elim_max = 3
[2025-10-16 16:02:16,314][write_input][INFO] n_subs = 3
[2025-10-16 16:02:16,314][write_input][INFO] subs_max = 3
[2025-10-16 16:02:16,314][write_input][INFO] target = random
[2025-10-16 16:02:16,314][write_input][INFO] end_point = (0.0, 0.0, 0.0)
[2025-10-16 16:02:16,314][write_input][INFO] show_max = 0.2
[2025-10-16 16:02:16,315][write_input][INFO] label_stable = True
[2025-10-16 16:02:16,315][write_input][INFO] vmax = 0.2
[2025-10-16 16:02:16,315][write_input][INFO] bottom_margin = 0.02
[2025-10-16 16:02:16,315][write_input][INFO] fig_format = svg
[2025-10-16 16:02:16,315][write_input][INFO]
[2025-10-16 16:02:16,315][write_input][INFO] [LAMMPS]
[2025-10-16 16:02:16,315][write_input][INFO] kpt_flag = False
[2025-10-16 16:02:16,315][write_input][INFO] force_gamma = False
[2025-10-16 16:02:16,315][write_input][INFO] lammps_infile = in.lmp
[2025-10-16 16:02:16,315][write_input][INFO] lammps_outfile = out.lmp
[2025-10-16 16:02:16,315][write_input][INFO] lammps_data = data.lmp
[2025-10-16 16:02:16,315][cryspy_init][INFO] # ---------- Initial structure generation
[2025-10-16 16:02:16,316][rs_gen][INFO] # ------ mindist
[2025-10-16 16:02:16,318][struc_util][INFO] Cu - Cu: 1.32
[2025-10-16 16:02:16,318][struc_util][INFO] Cu - Ag: 1.385
[2025-10-16 16:02:16,318][struc_util][INFO] Cu - Au: 1.34
[2025-10-16 16:02:16,318][struc_util][INFO] Ag - Ag: 1.45
[2025-10-16 16:02:16,318][struc_util][INFO] Ag - Au: 1.405
[2025-10-16 16:02:16,318][struc_util][INFO] Au - Au: 1.36
[2025-10-16 16:02:16,318][rs_gen][INFO] # ------ generate structures
[2025-10-16 16:02:16,395][gen_pyxtal][INFO] Structure ID      0: (4, 4, 0) Space group:  55 -->  69 Fmmm
[2025-10-16 16:02:16,448][gen_pyxtal][INFO] Structure ID      1: (2, 4, 3) Space group:  83 --> 123 P4/mmm
[2025-10-16 16:02:16,493][gen_pyxtal][INFO] Structure ID      2: (1, 4, 6) Space group:  25 -->  25 Pmm2
[2025-10-16 16:02:16,560][gen_pyxtal][INFO] Structure ID      3: (8, 2, 0) Space group: 109 --> 109 I4_1md
[2025-10-16 16:02:16,612][gen_pyxtal][INFO] Structure ID      4: (8, 7, 0) Space group:  22 -->  22 F222
[2025-10-16 16:02:16,650][gen_pyxtal][INFO] Structure ID      5: (2, 5, 5) Space group:   1 -->   1 P1
[2025-10-16 16:02:16,688][gen_pyxtal][INFO] Structure ID      6: (6, 2, 6) Space group:  84 --> 131 P4_2/mmc
[2025-10-16 16:02:16,727][gen_pyxtal][INFO] Structure ID      7: (1, 6, 6) Space group:  42 -->  42 Fmm2
[2025-10-16 16:02:16,738][gen_pyxtal][INFO] Structure ID      8: (0, 2, 6) Space group: 105 --> 105 P4_2mc
[2025-10-16 16:02:16,860][gen_pyxtal][INFO] Structure ID      9: (8, 7, 4) Space group: 111 --> 111 P-42m
[2025-10-16 16:02:16,901][gen_pyxtal][INFO] Structure ID     10: (3, 2, 7) Space group:  35 -->  35 Cmm2
[2025-10-16 16:02:17,007][gen_pyxtal][INFO] Structure ID     11: (6, 6, 4) Space group:  53 -->  53 Pmna
[2025-10-16 16:02:17,259][gen_pyxtal][INFO] Structure ID     12: (6, 6, 1) Space group: 200 --> 200 Pm-3
[2025-10-16 16:02:17,334][gen_pyxtal][INFO] Structure ID     13: (2, 1, 3) Space group:  79 --> 107 I4mm
[2025-10-16 16:02:17,359][gen_pyxtal][INFO] Structure ID     14: (8, 6, 2) Space group:  20 -->  20 C222_1
[2025-10-16 16:02:17,439][gen_pyxtal][INFO] Structure ID     15: (2, 8, 8) Space group:  21 -->  21 C222
[2025-10-16 16:02:17,506][gen_pyxtal][INFO] Structure ID     16: (5, 4, 3) Space group: 217 --> 217 I-43m
[2025-10-16 16:02:17,678][gen_pyxtal][INFO] Structure ID     17: (8, 7, 2) Space group: 164 --> 164 P-3m1
[2025-10-16 16:02:17,697][gen_pyxtal][INFO] Structure ID     18: (0, 3, 3) Space group: 154 --> 154 P3_221
[2025-10-16 16:02:17,710][gen_pyxtal][INFO] Structure ID     19: (2, 0, 8) Space group:  11 -->  11 P2_1/m
[2025-10-16 16:02:18,204][cryspy_init][INFO] Elapsed time for structure generation: 0:00:01.888822
[2025-10-16 16:02:18,208][ea_init][INFO] # ---------- Initialize evolutionary algorithm
[2025-10-16 16:02:18,208][ea_init][INFO] # ------ Generation 1
[2025-10-16 16:02:18,208][ea_init][INFO] 20 structures by random
[2025-10-16 16:02:19,428][cryspy_restart][INFO]
```
