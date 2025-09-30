NAME                           = "Freesurfer"
VERSION                        = 1
CONTRIBUTORS                   = ["ellis.langford.19@ucl.ac.uk"]
LAST_MOD_DATE                  = "29.09.2025"
VERBOSE                        = False

PARAMETERS = {
        "subject_id" : {
            "type"    : str,
            "help"    : "SubjectID of data to be analysed."
        },
        "input_im" : {
            "type"    : str,
            "help"    : "Path to a NIfTI image to be analysed."
        },
        "props_fpath" : {
            "type"       : str,
            "type_addit" : "json",
            "default"    : "",
            "help"       : "Path to a JSON-format plugin properties file (optional)."
        },
        "optional_flags" : {
            "type"       : str,
            "default"    : "",
            "help"       : "String of comma seperated flags to include in command. " +
                           "Default is none."
        },
        "parallelise" : {
            "type"    : bool,
            "default" : True,
            "help"    : "Run freesurfer processing in parallel. " +
                        "Default is True."
        },
        "no_randomness" : {
            "type"    : bool,
            "default" : True,
            "help"    : "Process with fixed random seeds to ensure reproducibility " +
                        "Default is True."
        },
        "big_vents" : {
            "type"    : bool,
            "default" : False,
            "help"    : "Use if subject has enlarged ventricles from atrophy. " +
                        "Default is False."
        },
        "large_FOV" : {
            "type"    : bool,
            "default" : False,
            "help"    : "Use if subject has a field of view > 256. " +
                        "Default is False."
        }
}