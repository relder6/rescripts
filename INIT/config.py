import os, sys, re

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
def get_flags():
    return {"USING_CURRENT_OFFSET": True,
            "USING_BOIL_CORR": True,
            "USING_CURRENT_CUT": True,
            "USING_DELTA_CORR": True,}


# ---------------------------------------------------------------------
# Data Cuts
# ---------------------------------------------------------------------
def get_data_cuts():
    return {"H_gtr_dp_min_cut": -8.0,
            "H_gtr_dp_max_cut": 8.0,
            "H_cer_npeSum_cut": 1.5,
            "H_cal_etottracknorm_cut": 0.7}

# ---------------------------------------------------------------------
# Data Cuts
# ---------------------------------------------------------------------
def get_common_values():
    return {"ebeam_4pass": 8.5831,
            "ebeam_5pass": 10.6716,
            "angle_4pass": 29.045,
            "angle_5pass": 16.75,
            "ebeam_3pass_phaseII": 6.449,
            "angle_3pass_phaseII": 31.278}

# ---------------------------------------------------------------------
# User Input Processing Logic
# ---------------------------------------------------------------------

# def get_common_run_inputs():
#     if len(sys.argv) == 4:
#         selected_run_type = sys.argv[1].strip().lower()
#         selected_beam_pass = sys.argv[2].strip()
#         selected_target = sys.argv[3].strip().lower()
#     else:
#         selected_run_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
#         if not selected_run_type:
#             selected_run_type = f"hmsdis"
#         selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
#         selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()


#     selected_beam_pass_to_energy_prefix = {"1": "2.","2": "4.","3": "6.","4": "8.","5": "10."}

#     beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)

#     if not beam_prefix:
#         print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
#         exit(1)

#     selected_target_shortcut_to_target_variable = {"al":"al","al13":"al","aluminum":"al",
#                                                    "c":"c","c12":"c","carbon":"c",
#                                                    "cu":"cu","cu29":"cu","copper":"cu",
#                                                    "opt1":"optics1","optics1":"optics1",
#                                                    "opt2":"optics2","optics2":"optics2",
#                                                    "d2":"ld2","ld2":"ld2",
#                                                    "h2":"lh2","lh2":"lh2",
#                                                    "hole":"hole","chole":"hole","c-hole":"hole",
#                                                    "dummy":"dummy","dum":"dummy",
#                                                    }

#     selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)

#     if not selected_target_shortname:
#         print(f"Unknown target: {selected_target}.  Please try again.")
#         exit(1)

#     selected_target_shortname_to_title_longname = {
#         "al":"Aluminum",
#         "c":"Carbon",
#         "cu":"Copper",
#         "opt1":"Optics1",
#         "opt2":"Optics2",
#         "ld2":"Deuterium",
#         "lh2":"Hydrogen",
#         "hole":"Carbon Hole",
#         "dummy":"Dummy"}

#     selected_target_titlename = selected_target_shortname_to_title_longname.get(selected_target_shortname)

#     selected_target_shortname_to_AZ = {"al":   (27, 13),
#                                        "c":    (12, 6),
#                                        "cu":   (64, 29),
#                                        "ld2":  (2,  1),
#                                        "lh2":  (1,  1),
#                                        "dummy": (0, 0)}
    
#     selected_target_A, selected_target_Z = selected_target_shortname_to_AZ.get(selected_target_shortname)

#     return selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z

# ---------------------------------------------------------------------
# User Input Processing Logic
# ---------------------------------------------------------------------
def parse_bins(nbins):
    if nbins is None or nbins.strip() == "":
        nbins = input("Indicate number of desired bins: ")
    nbins = int(nbins.strip())
    if not nbins:
        raise ValueError(f"Cannot evaluate number of bins: {nbins}")
    return nbins

def parse_run_type(selected_run_type):
    if selected_run_type is None or selected_run_type.strip() == "":
        selected_run_type = input("Enter desired run type (default HMSDIS): ").strip()
    selected_run_type = str(selected_run_type).strip().lower()
    if not selected_run_type:
        selected_run_type = "hmsdis"
    return selected_run_type

def parse_beam_pass(selected_beam_pass):
    selected_beam_pass_to_energy_prefix = {"1": "2.","2": "4.","3": "6.","4": "8.","5": "10."}
    if selected_beam_pass is None or selected_beam_pass.strip() == "":
        selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ")
    beam_prefix = selected_beam_pass_to_energy_prefix.get(str(selected_beam_pass).strip())
    if not beam_prefix:
        raise ValueError(f"Unknown pass: {selected_beam_pass}")
    return str(selected_beam_pass).strip(), beam_prefix


def parse_target(selected_target):
    if selected_target is None or selected_target.strip() == "":
        selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ")
    selected_target = str(selected_target).strip().lower()
    selected_target_shortcut_to_target_variable = {"al": "al", "aluminum": "al",
                                                   "c": "c", "carbon": "c",
                                                   "cu": "cu", "copper": "cu",
                                                   "ld2": "ld2", "d2": "ld2", "deuterium": "ld2",
                                                   "lh2": "lh2", "h2": "lh2", "hydrogen": "lh2",
                                                   "dummy": "dummy",
                                                   "dummy_up": "dummy_up", "up": "dummy_up",
                                                   "dummy_down": "dummy_down", "down": "dummy_down"}
    target_abbrev = selected_target_shortcut_to_target_variable.get(selected_target)
    if not target_abbrev:
        raise ValueError(f"Unknown target: {selected_target}")

    target_info = {"al": ("Aluminum", "Al", 27, 13),
                   "c": ("Carbon", "C", 12, 6),
                   "cu": ("Copper", "Cu", 64, 29),
                   "ld2": ("Deuterium", "LD2", 2, 1),
                   "lh2": ("Hydrogen", "LH2", 1, 1),
                   "dummy": ("Dummy", "Dummy", 0, 0),
                   "dummy_up": ("Dummy_Upstream", "Dummy_Up", 27, 13),
                   "dummy_down": ("Dummy_Downstream", "Dummy_Down", 27, 13)}
    target_longname, target_shortname, target_A, target_Z = target_info[target_abbrev]
    return target_abbrev, target_longname, target_shortname, target_A, target_Z

# def input_processing():
    
#     if len(sys.argv) == 4:
#         selected_run_type = sys.argv[1].strip().lower()
#         selected_beam_pass = sys.argv[2].strip()
#         selected_target = sys.argv[3].strip().lower()
        
#     else:
#         selected_run_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
#         selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
#         selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()
#         if not selected_run_type:
#             selected_run_type = f"hmsdis"

#     selected_beam_pass_to_energy_prefix = {"1": "2.","2": "4.","3": "6.","4": "8.","5": "10."}

#     beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)

#     if not beam_prefix:
#         print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
#         exit(1)
    
#     targets = {"c":   {"A": 12.0107,     "Z": 6,  "name": "Carbon",     "symbol": "C",   "aliases": ["c12","carbon"]},
#                "al":  {"A": 26.9815385,  "Z": 13, "name": "Aluminum",   "symbol": "Al",  "aliases": ["al13","aluminum"]},
#                "cu":  {"A": 63.546,      "Z": 29, "name": "Copper",     "symbol": "Cu",  "aliases": ["cu29","copper"]},
#                "fe":  {"A": 55.845,      "Z": 26, "name": "Iron",       "symbol": "Fe",  "aliases": ["iron"]},
#                "pb":  {"A": 207.2,       "Z": 82, "name": "Lead",       "symbol": "Pb",  "aliases": ["lead"]},
#                "ld2": {"A": 2.014101778, "Z": 1,  "name": "Deuterium",  "symbol": "D",   "aliases": ["d2","ld2"]},
#                "lh2": {"A": 1.008,       "Z": 1,  "name": "Hydrogen",   "symbol": "H",   "aliases": ["h2","lh2"]},}


