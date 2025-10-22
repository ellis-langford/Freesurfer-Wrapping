# Imports
import os
import sys
import subprocess
import nibabel as nib
import numpy as np
import ants
import trimesh
import shutil

class Core(object):
    """
    Technical Setup.
    """
    def __init__(self, plugin_obj):
        # Check all expected attributed are present
        to_inherit = ["config", "utils", "helpers", "parameters",
                      "base_dir", "fs_outputs", "output_dir", 
                      "image_dir", "log_dir", "interim_dir", 
                      "subject_id", "input_im", 
                      "freesurfer_outputs", "freesurfer_command",
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
        segmentation = os.path.join(self.fs_outputs, "mri", "aseg.mgz")
        output_im    = os.path.join(self.fs_outputs, "mri", "T1.mgz")
        if not os.path.exists(segmentation):
            self.helpers.errors(f"Freesurfer has not produced a segmentation at {segmentation}" +
                                f"- please check log file at {self.freesurfer_log}")
        elif not os.path.exists(output_im):
            self.helpers.errors(f"Freesurfer has not produced an output image at {output_im} - " +
                                f"please check log file at {self.freesurfer_log}")
        else:
            self.helpers.plugin_log("Freesurfer run successfully")
            
        return

    def convert_T1(self):
        """
        Convert freesurfer T1.mgz to .nii.gz
        """
        self.helpers.verbose_log("Converting T1.mgz to T1.nii.gz")

        # Define image and log paths
        T1 = os.path.join(self.fs_outputs, "mri", "T1.mgz")
        self.T1_out  = os.path.join(self.fs_outputs, "mri", "T1.nii.gz")
        conversion_log = os.path.join(self.log_dir, "T1_conversion.log")

        # Convert
        with open(conversion_log, "w") as outfile:
            subprocess.run(["bash", "-c",
                            self.freesurfer_source + "mri_convert " +
                            f"{T1} {self.T1_out}"],
                            stdout=outfile,
                            stderr=subprocess.STDOUT,
                            env=self.freesurfer_env)
            
        # Check if conversion successful
        if not os.path.exists(self.T1_out):
            self.helpers.errors("Conversion of T1.mgz to T1.nii.gz failed - " +
                               f"please check log file at {conversion_log}")
        else:
            shutil.copy(self.T1_out, os.path.join(self.output_dir, "T1.nii.gz"))
            self.helpers.verbose_log("Conversion of T1.mgz to T1.nii.gz successful")
        
        return

    def convert_seg(self):
        """
        Convert freesurfer aseg.mgz to .nii.gz
        """
        self.helpers.verbose_log("Converting aseg.mgz to aseg.nii.gz")

        # Define image and log paths
        seg = os.path.join(self.fs_outputs, "mri", "aseg.mgz")
        self.seg_out  = os.path.join(self.fs_outputs, "mri", "aseg.nii.gz")
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
            self.helpers.verbose_log("Conversion of aseg.mgz to aseg.nii.gz successful")
        
        return

    def binarise(self, region):
        """
        Binarise to extract required freesurfer labels
        """
        self.helpers.verbose_log(f"Binarising {region} segmentation")

        # Define image and log paths
        subcortical_seg = os.path.join(self.fs_outputs, "mri", "aseg.nii.gz")
        bin_out  = os.path.join(self.interim_dir, f"{region}", f"{region}_bin.nii.gz")
        os.makedirs(os.path.join(self.interim_dir, region), exist_ok=True)
        binarise_log = os.path.join(self.log_dir, f"binarise.log")

        labels = {"ventricles"     : "24 4 5 14 15 43 44 213",
                  "cerebellum_L"   : "6 7 8",
                  "cerebellum_R"   : "45 46 47",
                  "cerebellumWM_L" : "7",
                  "cerebellumWM_R" : "46",
                  "brainstem"      : "16 170 171 172 173 174 175 177 178 179 71000 71010",
                  "cerebrum_L"     : "2 3 10 11 12 13 17 18 19 20 26 28",
                  "cerebrum_R"     : "41 42 49 50 51 52 53 54 55 56 58 60",
                  "cerebrumWM_L"   : "2 78",
                  "cerebrumWM_R"   : "41 79",
                  "wholebrain"     : "6 7 8 16 45 46 47 192 "
                                     "24 4 5 14 15 43 44 213 "
                                     "2 3 10 11 12 13 17 18 19 20 26 28 "
                                     "41 42 49 50 51 52 53 54 55 56 58 60"
                 }

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
        atlas_t1 = "/app/atlas/mni_icbm152_atlas_t1.nii.gz"
        atlas_labels = "/app/atlas/mni_icbm152_atlas_t1.nii.gz"
        atlas_labels_out = os.path.join(self.interim_dir, "mni_icbm152_labels_subjectspace.nii.gz")
        brainstem_seg = os.path.join(self.interim_dir, "brainstem", "brainstem_bin.nii.gz")
        input_image = self.input_im if self.input_im else os.path.join(self.fs_outputs, "mri", "T1.nii.gz")
        
        # Register atlas T1 to subject T1 space
        registration = ants.registration(fixed=ants.image_read(input_image), 
                                         moving=ants.image_read(atlas_t1), 
                                         type_of_transform="Affine")

        # Apply transform to atlas labels
        transformed_labels = ants.apply_transforms(fixed=ants.image_read(input_image),
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
        input_im = os.path.join(self.interim_dir, f"{region}", "smooth")

        # Define output and log paths
        geo_out = os.path.join(self.interim_dir, f"{region}", f"{region}.stl")
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

    def clean_stl(self, region):
        """
        Fix issues in stl file (e.g. non-normal orientations)
        """
        self.helpers.verbose_log(f"Cleaning {region} .stl")
        input_im = os.path.join(self.interim_dir, f"{region}", f"{region}.stl")

        # Define output and log paths
        geo_out = os.path.join(self.output_dir, f"{region}.stl")

        # Clean
        mesh = trimesh.load(input_im) # Load mesh
        components = mesh.split(only_watertight=False) # Split into connected components
        largest = max(components, key=lambda c: c.area) # Keep largest shell
        
        # Remove small disconnected patches
        cleaned_components = [
            c for c in components if c.area > 0.05 * largest.area
        ]
        
        mesh_clean = trimesh.util.concatenate(cleaned_components)
        
        # Fill holes
        if not mesh_clean.is_watertight:
            mesh_clean.fill_holes()
        
        # Export cleaned STL
        mesh_clean.export(geo_out)
            
        # Check if conversion successful
        if not os.path.exists(geo_out):
            self.helpers.errors(f"Cleaning of {region} .stl failed")
        else:
            self.helpers.verbose_log(f"Cleaning of {region} .stl successful")
        
        return

    def prepare_geometry(self, region):
        """
        Binarises, tessellates, smooths, converts and cleans volume to .stl
        """
        self.helpers.verbose_log(f"Preparing {region} geometry")

        # Process geometry
        self.tessellate(region)
        self.smooth(region)
        self.convert_to_stl(region)
        self.clean_stl(region)

        # Define expected outputs
        surf_out = os.path.join(self.interim_dir, f"{region}", "surf")
        smooth_out = os.path.join(self.interim_dir, f"{region}", "smooth")
        conv_out = os.path.join(self.interim_dir, f"{region}", f"{region}.stl")
        geo_out = os.path.join(self.output_dir, f"{region}.stl")
            
        # Check if geometry generation successful
        outputs = [surf_out, smooth_out, conv_out, geo_out]
        for output in outputs:
            if not os.path.exists(output):
                self.helpers.errors(f"Geometry generation failed - " +
                                    f"required output missing ({output})")

        self.helpers.verbose_log(f"Generation of {region} geometry successful")
        
        return

    def generate_global_surface(self):
        """
        Generate global mesh file
        (all regions minus ventricles)
        """
        # Load global binary and ventricle mask
        wholebrain_seg = nib.load(os.path.join(self.interim_dir, "wholebrain", "wholebrain_bin.nii.gz"))
        vent_seg = nib.load(os.path.join(self.interim_dir, "ventricles", "ventricles_bin.nii.gz"))
        wholebrain_data = wholebrain_seg.get_fdata()
        vent_data = vent_seg.get_fdata()
    
        # Subtract ventricles (make sure masks are binary)
        result_data = np.where((wholebrain_data > 0) & (vent_data == 0), 1, 0)
    
        # Step 4: Save result
        result_img = nib.Nifti1Image(result_data.astype(np.uint8), affine=wholebrain_seg.affine)
        os.makedirs(os.path.join(self.interim_dir, "global"), exist_ok=True)
        result_out_fpath = os.path.join(self.interim_dir, "global", "global_bin.nii.gz")
        nib.save(result_img, result_out_fpath)


    def run(self):
        """
        Begin core processing.
        """
        self.helpers.plugin_log(f"Starting core processing at {self.helpers.now_time()}")

        # Record parameters
        self.helpers.log_options(self.parameters)
        regions = self.parameters["regions"].split(",")

        # Run Freesurfer
        if self.input_im:
            self.helpers.plugin_log("Running Freesurfer")
            self.run_freesurfer()

        if not self.parameters["segmentations"]:
            # Convert aseg file to NIfTI
            self.helpers.plugin_log("Converting aseg and T1 file to NIfTI")
            self.convert_seg()
            self.convert_T1()

            # Create region binary files
            self.helpers.plugin_log("Creating region binary files")
            
            # Join brainstem label - split in later steps
            if "brainstem_L" in regions or "brainstem_R" in regions:
                binarise_regions = [r for r in regions if r not in ("brainstem_L", "brainstem_R")]
                binarise_regions.append("brainstem")

            # Process regions
            for region in binarise_regions:
                self.binarise(region)
    
            # Register MNI-ICBM152 atlas labels to subject space
            self.helpers.plugin_log("Registering atlas labels to subject space")
            self.register_mni_atlas()
    
            # Split the brainstem into L&R
            self.helpers.plugin_log("Splitting brainstem into L&R")
            self.split_brainstem()

        else:
            for region in regions:
                self.utils.copy(os.path.join(self.image_dir, f"{region}_bin.nii.gz"),
                                os.path.join(self.interim_dir, region, f"{region}_bin.nii.gz"))

        # Create other ROI geometries
        self.helpers.plugin_log(f"Creating region surface files")
        for region in regions:
            if self.parameters["segmentations"]:
                self.prepare_geometry(region)
            else:
                 if region not in ["wholebrain", "ventricles"]:
                     self.prepare_geometry(region)

        # If segmentation mode with both wholebrain and ventricles in region
        # list, create global segmentation (all minus ventricles)
        if self.parameters["segmentations"]:
            if "wholebrain" in regions and "ventricles" in regions:
                self.generate_global_surface()
                self.prepare_geometry("global")

        self.helpers.plugin_log(f"Core analysis completed at {self.helpers.now_time()}")
        self.helpers.log_success()

        return