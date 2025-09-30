# Imports
import os
import sys
import subprocess
import nibabel as nib
import numpy as np
import ants

class Core(object):
    """
    Technical Setup.
    """
    def __init__(self, plugin_obj):
        # Check all expected attributed are present
        to_inherit = ["config", "utils", "helpers", "parameters",
                      "base_dir", "output_dir", "image_dir", "log_dir",
                      "interim_dir", "geo_dir", "subject_id", 
                      "input_im", "freesurfer_command",
                      "freesurfer_source", "freesurfer_env"]
        for attr in to_inherit:
            try:
                setattr(self, attr, getattr(plugin_obj, attr))
            except AttributeError as e:
                print(f"Attribute Error - {e}")
                sys.exit(1)

    def run_freesurfer(self):
        """
        Run Freesurfer.
        """
        self.helpers.plugin_log("Starting execution")
        self.freesurfer_log     = os.path.join(self.log_dir, "freesurfer.log")
        self.freesurfer_command = self.freesurfer_source + self.freesurfer_command
        self.helpers.plugin_log(f"Freesurfer command: {self.freesurfer_command}")
        with open(self.freesurfer_log, "w") as outfile:
            freesurfer_sub = subprocess.run(["bash", "-c",
                                             self.freesurfer_command],
                                             stdout=outfile,
                                             stderr=subprocess.STDOUT,
                                             env=self.freesurfer_env)
        if freesurfer_sub.returncode != 0:
            self.helpers.errors(f"Freesurfer execution returned non-zero exit status - " +
                                f"please check log file at {self.freesurfer_log}")

        # Check required outputs have been produced
        segmentation = os.path.join(self.output_dir, "mri", "aseg.mgz")
        output_im    = os.path.join(self.output_dir, "mri", "T1.mgz")
        if not os.path.exists(segmentation):
            self.helpers.errors(f"Freesurfer has not produced a segmentation at {segmentation}" +
                                f"- please check log file at {self.freesurfer_log}")
        elif not os.path.exists(output_im):
            self.helpers.errors(f"Freesurfer has not produced an output image at {output_im} - " +
                                f"please check log file at {self.freesurfer_log}")
        else:
            self.helpers.plugin_log("Freesurfer run successfully")
            
        return

    def convert_seg(self):
        """
        Convert freesurfer aseg.mgz to .nii.gz
        """
        self.helpers.verbose_log("Converting aseg.mgz to aseg.nii.gz")

        # Define image and log paths
        seg = os.path.join(self.output_dir, "mri", "aseg.mgz")
        self.seg_out  = os.path.join(self.output_dir, "mri", "aseg.nii.gz")
        conversion_log = os.path.join(self.log_dir, "aseg_conversion.log")

        # Convert
        with open(conversion_log, "w") as outfile:
            subprocess.run(["bash", "-c",
                            self.freesurfer_source + "mri_convert " +
                            f"{seg} {self.seg_out}"],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                            env=self.freesurfer_env)
            
        # Check if conversion successful
        if not os.path.exists(self.seg_out):
            self.helpers.errors("Conversion of aseg.mgz to aseg.nii.gz failed - " +
                               f"please check log file at {conversion_log}")
        else:
            self.helpers.verbose_log("Conversion of aseg.mgz to nii successful")
        
        return

    def binarise(self, region):
        """
        Binarise to extract required freesurfer labels
        """
        self.helpers.verbose_log(f"Binarising {region} segmentation")

        # Define image and log paths
        subcortical_seg = os.path.join(self.output_dir, "mri", "aseg.nii.gz")
        bin_out  = os.path.join(self.interim_dir, f"{region}", f"{region}_bin.nii.gz")
        os.makedirs(os.path.join(self.interim_dir, region), exist_ok=True)
        binarise_log = os.path.join(self.log_dir, f"binarise.log")

        labels = {"ventricles"     : "24 4 5 14 15 43 44 213",
                  "cerebellum_L"   : "8",
                  "cerebellum_R"   : "47",
                  "cerebellumWM_L" : "7",
                  "cerebellumWM_R" : "46",
                  "brainstem"      : "16 170 171 172 173 174 175 177 178 179 71000 71010"}

        # Convert
        with open(binarise_log, "w") as outfile:
            subprocess.run(["bash", "-c",
                            self.freesurfer_source + "mri_binarize " +
                            f"--i {subcortical_seg} --match {labels[region]}" +
                            f"--inv --o {bin_out}"],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                            env=self.freesurfer_env)
            
        # Check if conversion successful
        if not os.path.exists(bin_out):
            self.helpers.errors(f"Binarisation of {region} segmentation failed - " +
                                f"please check log file at {binarise_log}")
        else:
            self.helpers.verbose_log(f"Binarisation of {region} segmentation successful")
        
        return

    def register_mni_atlas(self):
        """
        Register MNI-ICBM152 CerebrA atlas labels to subject space.
    
        Reference:
            Manera AL, Dadar M, Fonov V, Collins DL. (2020).
            CerebrA, registration and manual label correction of Mindboggle-101 atlas
            for MNI-ICBM152 template. Scientific Data, 7, 237.
            https://doi.org/10.1038/s41597-020-00564-0
        """
        # Set paths
        atlas_t1 = "/app/atlas/mni_icbm152_atlas_t1.nii"
        atlas_labels = "/app/atlas/mni_icbm152_atlas_t1.nii"
        atlas_labels_out = os.path.join(self.interim_dir, "mni_icbm152_labels_subjectspace.nii")
        brainstem_seg = os.path.join(self.interim_dir, "brainstem", "brainstem_bin.nii.gz")
        
        # Register atlas T1 to subject T1 space
        registration = ants.registration(fixed=ants.image_read(self.input_im), 
                                         moving=ants.image_read(atlas_t1), 
                                         type_of_transform="Affine")

        # Apply transform to atlas labels
        transformed_labels = ants.apply_transforms(fixed=ants.image_read(self.input_im),
                                                   moving=ants.image_read(atlas_labels),
                                                   transformlist=[registration["fwdtransforms"][0]],
                                                   interpolator="nearestNeighbor")

        # Resample labels so dimensions match brainstem seg
        resampled = ants.resample_image_to_target(image=transformed_labels,
                                                  target=ants.image_read(brainstem_seg),
                                                  interp_type="nearestNeighbor")
                
        # Save transformed labels
        ants.image_write(resampled, atlas_labels_out)
        self.atlas_in_subj = atlas_labels_out

        # Check outputs
        if not os.path.exists(atlas_labels_out):
            self.helpers.errors("Transformation of atlas labels to subject space failed")
        else:
            self.helpers.verbose_log("Transformation of atlas labels to subject space successful")

        return

    def split_brainstem(self):
        """
        Split the binarised brainstem into left and right hemispheres
        """
        self.helpers.verbose_log("Splitting brainstem labels into left/right")
        
        # Load brainstem segmentation
        brainstem_seg = nib.load(os.path.join(self.interim_dir, "brainstem", "brainstem_bin.nii.gz"))
        brainstem_data = brainstem_seg.get_fdata().astype(np.uint8)
        affine = brainstem_seg.affine
        
        # Load atlas labels
        atlas_img = nib.load(self.atlas_in_subj)
        atlas_data = atlas_img.get_fdata().astype(int)
        
        # Create masks
        mask = brainstem_data > 0
        left_mask = np.zeros_like(brainstem_data, dtype=np.uint8)
        right_mask = np.zeros_like(brainstem_data, dtype=np.uint8)
        
        # Assign voxels that overlap atlas left/right IDs
        left_vox = np.isin(atlas_data, [62]) & mask
        right_vox = np.isin(atlas_data, [11]) & mask
        left_mask[left_vox] = 1
        right_mask[right_vox] = 1
        
        # Find leftover voxels (in aseg but not in atlas L/R)
        assigned = left_vox | right_vox
        leftovers = mask & (~assigned)
        
        if np.any(leftovers):
            coords = np.array(np.nonzero(leftovers)).T
            ras_coords = nib.affines.apply_affine(affine, coords)
        
            for (i, j, k), ras in zip(coords, ras_coords):
                if ras[0] < 0:   # Left side
                    left_mask[i, j, k] = 1
                else:            # Right side
                    right_mask[i, j, k] = 1
        
        # Save outputs
        left_output = os.path.join(self.interim_dir, "brainstem_L", "brainstem_L_bin.nii.gz")
        right_output = os.path.join(self.interim_dir, "brainstem_R", "brainstem_R_bin.nii.gz")
        for file in [left_output, right_output]:
            _dir = os.path.dirname(file)
            os.makedirs(_dir, exist_ok=True)
            
        nib.save(nib.Nifti1Image(left_mask, affine, brainstem_seg.header), left_output)
        nib.save(nib.Nifti1Image(right_mask, affine, brainstem_seg.header), right_output)

        # Check outputs
        if not os.path.exists(left_output):
            self.helpers.errors(f"Splitting of brainstem region failed")
        elif not os.path.exists(right_output):
            self.helpers.errors(f"Splitting of brainstem region failed")
        else:
            self.helpers.verbose_log(f"Splitting of brainstem successful")

        return

    def tessellate(self, region):
        """
        Tessellate to create surface from input volume
        """
        self.helpers.verbose_log(f"Tessellating {region}")

        # Define image and log paths
        bin_data = os.path.join(self.interim_dir, f"{region}", f"{region}_bin.nii.gz")
        surf_out = os.path.join(self.interim_dir, f"{region}", "surf")
        tessellate_log = os.path.join(self.log_dir, f"tessellate_{region}.log")

        # Convert
        with open(tessellate_log, "w") as outfile:
            subprocess.run(["bash", "-c",
                            self.freesurfer_source + "mri_tessellate " +
                            f"{bin_data} 1 {surf_out}"],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                            env=self.freesurfer_env)
            
        # Check if conversion successful
        if not os.path.exists(surf_out):
            self.helpers.errors(f"Tessellation of {region} failed - " +
                                f"please check log file at {tessellate_log}")
        else:
            self.helpers.verbose_log(f"Tessellation of {region} successful")
        
        return

    def smooth(self, region):
        """
        Smooths the tessellation of region surface
        """
        self.helpers.verbose_log(f"Smoothing {region}")

        # Define image and log paths
        surfs = os.path.join(self.interim_dir, f"{region}", "surf")
        smooth_out = os.path.join(self.interim_dir, f"{region}", "smooth")
        smooth_log = os.path.join(self.log_dir, f"smooth_{region}.log")

        # Convert
        with open(smooth_log, "w") as outfile:
            subprocess.run(["bash", "-c",
                            self.freesurfer_source + "mris_smooth " +
                            f"{surfs} {smooth_out}"],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                            env=self.freesurfer_env)
            
        # Check if conversion successful
        if not os.path.exists(smooth_out):
            self.helpers.errors(f"Smoothing of {region} failed - " +
                                f"please check log file at {smooth_log}")
        else:
            self.helpers.verbose_log(f"Smoothing of {region} successful")
        
        return

    def convert_to_stl(self, region):
        """
        Converts geometry file to .stl
        """
        self.helpers.verbose_log(f"Converting {region} to .stl")
        paths = {
            "cerebrum_L" : os.path.join(self.output_dir, "surf", "lh.pial"),
            "cerebrum_R" : os.path.join(self.output_dir, "surf", "rh.pial"),
            "cerebrumWM_L" : os.path.join(self.output_dir, "surf", "lh.smoothwm"),
            "cerebrumWM_R" : os.path.join(self.output_dir, "surf", "rh.smoothwm")
        }

        # Define input images
        if region in paths:
            input_im = paths[region]
        else:
            input_im = os.path.join(self.interim_dir, f"{region}", "smooth")

        # Define output and log paths
        geo_out = os.path.join(self.geo_dir, f"{region}.stl")
        conversion_log = os.path.join(self.log_dir, f"{region}_conversion.log")

        # Convert
        with open(conversion_log, "w") as outfile:
            subprocess.run(["bash", "-c",
                            self.freesurfer_source + "mris_convert " +
                            f"{input_im} {geo_out}"],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                            env=self.freesurfer_env)
            
        # Check if conversion successful
        if not os.path.exists(geo_out):
            self.helpers.errors(f"Conversion {region} to .stl failed - " +
                                f"please check log file at {conversion_log}")
        else:
            self.helpers.verbose_log(f"Conversion of {region} to .stl successful")
        
        return

    def prepare_geometry(self, region):
        """
        Binarises, tessellates, smooths and converts volume to .stl
        """
        self.helpers.verbose_log(f"Preparing {region} geometry")

        # Process geometry
        self.tessellate(region)
        self.smooth(region)
        self.convert_to_stl(region)

        # Define expected outputs
        surf_out = os.path.join(self.interim_dir, f"{region}", "surf")
        smooth_out = os.path.join(self.interim_dir, f"{region}", "smooth")
        geo_out = os.path.join(self.geo_dir, f"{region}.stl")
            
        # Check if geometry generation successful
        outputs = [surf_out, smooth_out, geo_out]
        for output in outputs:
            if not os.path.exists(output):
                self.helpers.errors(f"Geometry generation failed - " +
                                    f"required output missing ({output})")

        self.helpers.verbose_log(f"Generation of {region} geometry successful")
        
        return

    def run(self):
        """
        Begin core processing.
        """
        self.helpers.plugin_log(f"Starting core processing at {self.helpers.now_time()}")

        # Record parameters
        self.helpers.log_options(self.parameters)

        # Run Freesurfer
        self.helpers.plugin_log("Running Freesurfer")
        self.run_freesurfer()

        # Convert aseg file to NIfTI
        self.helpers.plugin_log("Converting aseg file to NIfTI")
        self.convert_seg()

        # Create cerebrum geometry
        self.helpers.plugin_log("Creating cerebrum surface files")
        regions = ["cerebrum_L", "cerebrum_R", "cerebrumWM_L", "cerebrumWM_R"]
        for region in regions:
            self.convert_to_stl(region)

        # Create region binary files
        self.helpers.plugin_log(f"Creating region binary files")
        regions = ["brainstem", "cerebellum_L", "cerebellum_L",
                   "cerebellumWM_L", "cerebellumWM_R", "ventricles"]
        for region in regions:
            self.binarise(region)

        # Register MNI-ICBM152 atlas labels to subject space
        self.helpers.plugin_log("Registering atlas labels to subject space")
        self.register_mni_atlas()

        # Split the brainstem into L&R
        self.helpers.plugin_log("Splitting brainstem into L&R")
        self.split_brainstem()

        # Create other ROI geometries
        self.helpers.plugin_log(f"Creating region surface files")
        regions = ["brainstem_L", "brainstem_R", "cerebellum_L", "cerebellum_L",
                   "cerebellumWM_L", "cerebellumWM_R", "ventricles"]
        for region in regions:
            self.prepare_geometry(region)

        self.helpers.plugin_log(f"Core analysis completed at {self.helpers.now_time()}")

        return