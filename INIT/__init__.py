import os, sys, re

# =====================================================================
# Handling user inputs
# =====================================================================

def get_common_run_inputs():
    if len(sys.argv) == 4:
        selected_run_type = sys.argv[1].strip().lower()
        selected_beam_pass = sys.argv[2].strip()
        selected_target = sys.argv[3].strip().lower()
    else:
        selected_run_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
        if not selected_run_type:
            selected_run_type = f"hmsdis"
        selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
        selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()


    selected_beam_pass_to_energy_prefix = {"1": "2.","2": "4.","3": "6.","4": "8.","5": "10."}

    beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)

    if not beam_prefix:
        print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
        exit(1)

    selected_target_shortcut_to_target_variable = {"al":"al","al13":"al","aluminum":"al",
                                                   "c":"c","c12":"c","carbon":"c",
                                                   "cu":"cu","cu29":"cu","copper":"cu",
                                                   "opt1":"optics1","optics1":"optics1",
                                                   "opt2":"optics2","optics2":"optics2",
                                                   "d2":"ld2","ld2":"ld2",
                                                   "h2":"lh2","lh2":"lh2",
                                                   "hole":"hole","chole":"hole","c-hole":"hole",
                                                   "dummy":"dummy","dum":"dummy",
                                                   }

    selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)

    if not selected_target_shortname:
        print(f"Unknown target: {selected_target}.  Please try again.")
        exit(1)

    return selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname
