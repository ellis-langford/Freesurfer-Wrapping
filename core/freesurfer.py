## Top-level wrapper
# Imports
import os
import sys
import shutil
from core import Core

# Add the utils directory to the Python path
utils_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.append(utils_path)

# Import custom modules
try:
    from base_cog import BaseCog
    from utils import Utils
    from helpers import Helpers

except Exception as e:
    print(f"Failed to import toolkit modules from /app/utils - {e}")
    sys.exit(1)

# Create a class for processing
class Freesurfer(BaseCog):
    def __init__(self, **kwargs):
        """
        Instantiate the Freesurfer class.
        """
        super().__init__(**kwargs)
        
        # Instantiate custom modules
        self.utils = Utils()
        self.helpers = Helpers()

        # Load parameters from CLI or properties file
        self.load_parameters()

    def copy_inputs(self):
        """
        Check inputs and copy to working directory.
        """
        input_im = self.get_parameter("input_im")
        freesurfer_outputs = self.get_parameter("freesurfer_outputs")
        segmentations = self.get_parameter("segmentations")
        regions = self.get_parameter("regions")
        
        # Check for valid inputs and copy
        # Run full FreeSurfer pipeline
        if input_im:
            if not os.path.isfile(input_im):
                self.helpers.errors(f"No valid image found at {input_im}")
            else:
                self.utils.copy(input_im, os.path.join(self.image_dir, "image.nii.gz"))
        # Post processing only
        elif freesurfer_outputs:
            if not os.path.isdir(freesurfer_outputs):
                self.helpers.errors(f"No valid Freesurfer output directory found at {freesurfer_outputs}")
            else:
                shutil.copytree(freesurfer_outputs, os.path.join(self.base_dir, "fs_outputs", self.subject_id))
        # Post processing of segmentations only
        elif segmentations:
            segmentation_paths = segmentations.split(",")
            regions = regions.split(",")
            for region in regions:
                path = next((p for p in segmentation_paths if region in os.path.basename(p)), None)
                if not path:
                    self.helpers.errors(f"No matching path in segmentations found for {region}")
                elif not os.path.isfile(path):
                    self.helpers.errors(f"No valid segmentation found at {path}")
                else:
                    self.utils.copy(path, os.path.join(self.image_dir, f"{region}_bin.nii.gz"))
        else:
            self.helpers.errors(f"A T1, directory containing Freesurfer outputs or "
                                f"binary segmentation .nii.gz image must be provided")
            
        # Assign variables
        self.input_im = os.path.join(self.image_dir, "image.nii.gz") if input_im else None
        self.freesurfer_outputs = os.path.join(self.base_dir, "fs_outputs", self.subject_id) if freesurfer_outputs else None

        return

    def build_freesurfer_command(self):
        """
        Build freesurfer command.
        """
        # Define command
        self.freesurfer_source  = "source $FREESURFER_HOME/SetUpFreeSurfer.sh && "
        self.freesurfer_command = f"recon-all -i {self.input_im} -subject {self.subject_id} -all"

        # Add extra flags
        for _input, tag in {"parallelise"        : "parallel",
                            "no_randomness"      : "norandomness",
                            "big_vents"          : "bigventricles",
                            "large_FOV"          : "cw256"}.items():

            if self.get_parameter(_input):
                self.freesurfer_command += f" -{tag}"

        # Add optional custom flags
        if self.get_parameter("optional_flags"):
            optional_flags_str = self.get_parameter("optional_flags")
            optional_flags = optional_flags_str.split(",")
            for flag in optional_flags:
                self.freesurfer_command += f" -{flag}"

        return

    def core(self):
        """
        Core processing
        """
        self.helpers.plugin_log(f"Starting {self.config['NAME']} at {self.helpers.now_time()}")
        
        # Tidy up log files from previous runs
        self.helpers.tidy_up_logs()

        # Define subjectID
        self.subject_id = self.get_parameter("subject_id")

        # Directories
        self.image_dir     = os.path.join(self.base_dir, "inputs")
        self.interim_dir   = os.path.join(self.base_dir, "interim")
        self.output_dir       = os.path.join(self.base_dir, "outputs")
        self.fs_output_dir = os.path.join(self.base_dir, "fs_outputs")
        self.fs_outputs    = os.path.join(self.base_dir, "fs_outputs", self.subject_id)
        self.log_dir   = os.path.join(self.base_dir, "logs")

        for _dir in [self.image_dir, self.interim_dir, self.output_dir, 
                     self.fs_output_dir, self.log_dir]:
            shutil.rmtree(_dir, ignore_errors=True)
            os.makedirs(_dir, exist_ok=True)
        
        # Record parameters
        self.helpers.log_options(self.parameters)
        
        # Copy inputs
        self.helpers.plugin_log("Copying inputs")
        self.copy_inputs()

        # Initialise Freesurfer variables and build command
        self.freesurfer_env                 = os.environ.copy()
        self.freesurfer_env["SUBJECTS_DIR"] = self.fs_output_dir
        self.build_freesurfer_command()

        # Core processing
        proc = Core(self)
        proc.run()

        # Complete
        self.helpers.plugin_log(f"{self.config['NAME']} completed at {self.helpers.now_time()}")
        self.helpers.log_success()

# Set off pipeline object
if __name__ == "__main__":
    # Init
    plugin = Freesurfer()
    # Execute
    plugin.core()