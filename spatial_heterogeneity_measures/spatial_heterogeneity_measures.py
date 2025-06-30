#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file only contains the code for functions used to derive quantitative spatial heterigeneity measures.

@author: P-K-W
"""

import nibabel as nib
import numpy as np
import math
from scipy.ndimage import convolve

#%% [markdown]
# ##These are the functiosn that were used to calculate the spatial heterogeneity measures
#

#%%

#Load nifti file and returns its data as a Numpy array
def load_nifti(nifti_path):
    return nib.load(nifti_path).get_fdata()

#Converts z-score mmap into a  binary mask
def binarize_zmap(zmap, threshold=3.09):
    return zmap > threshold  # Returns boolean mask

#GOF calculation
def compute_gof(mask1, mask2):

    # Initialize list for GoF scores
    gof_scores = []

    # Loop through each mode (4th dimension)
    for mode in range(mask1.shape[3]):  
        mode_data = mask1[:, :, :, mode]  # Extract the mode

        # Compute overlapping activated voxels
        overlapped_voxels = np.logical_and(mode_data, mask2)

        # Compute non-overlapping activated voxels (activated in mode but NOT in SN_mask)
        non_overlapped_voxels = np.logical_and(mode_data, np.logical_not(mask2))

        # Compute total activated voxels in the current mode
        total_activated_voxels = np.sum(mode_data)

        # Compute GoF score
        gof_score = (np.sum(overlapped_voxels) - np.sum(non_overlapped_voxels)) / total_activated_voxels
        
        # Store result
        gof_scores.append(gof_score)

    return gof_scores

#Median GOF calculation
def  median_gof(zmap, mask2, threshold=3.09):
    # Initialize list for GoF scores
    median_gof_scores = []

    # Loop through each mode (4th dimension)
    for mode in range(zmap.shape[3]):  
        mode_data = zmap[..., mode]  # Extract the 3D z score map for the current mode
        
        # Identify activated voxels (z > threshold)
        activated = mode_data > threshold
        
        # Compute overlapping activated voxels (activated and within mask2)
        overlapped_voxels = activated & mask2
        
        # Compute non-overlapping activated voxels (activated but outside mask2)
        non_overlapped_voxels = activated & (~mask2)
        
        # Compute the average z score for voxels overlapping mask2
        if np.sum(overlapped_voxels) > 0:
            med_overlap = np.median(mode_data[overlapped_voxels])
        else:
            med_overlap = np.nan  # Alternatively, use 0 or handle as needed
        
        # Compute the average z score for voxels outside mask2
        if np.sum(non_overlapped_voxels) > 0:
            med_nonoverlap = np.median(mode_data[non_overlapped_voxels])
        else:
            med_nonoverlap = np.nan  # Alternatively, use 0 or handle as needed
        
        # Compute the difference
        score = med_overlap - med_nonoverlap
        
        # Store the result for this mode
        median_gof_scores.append(score)
    
    return median_gof_scores

#Weighted GOF calculation
def  weighted_gof(zmap, mask2, threshold=3.09):
    # Initialize list for GoF scores
    weighted_gof_scores = []

    # Loop through each mode (4th dimension)
    for mode in range(zmap.shape[3]):  
        mode_data = zmap[..., mode]  # Extract the 3D z score map for the current mode
        
        # Identify activated voxels (z > threshold)
        activated = mode_data > threshold
        
        # Compute overlapping activated voxels (activated and within mask2)
        overlapped_voxels = activated & mask2
        
        # Compute non-overlapping activated voxels (activated but outside mask2)
        non_overlapped_voxels = activated & (~mask2)
        
        # Compute the average z score for voxels overlapping mask2
        if np.sum(overlapped_voxels) > 0:
            med_overlap = np.median(mode_data[overlapped_voxels])
        else:
            med_overlap = np.nan  # Alternatively, use 0 or handle as needed
        
        # Compute the average z score for voxels outside mask2
        if np.sum(non_overlapped_voxels) > 0:
            med_nonoverlap = np.median(mode_data[non_overlapped_voxels])
        else:
            med_nonoverlap = np.nan  # Alternatively, use 0 or handle as needed
        
        # Compute total activated voxels in the current mode
        total_activated_voxels = np.sum(activated)
        
        # Compute the weighted gof score
        score = med_overlap * (np.sum(overlapped_voxels) / total_activated_voxels) - med_nonoverlap * (np.sum(non_overlapped_voxels) / total_activated_voxels)

        # Store the result for this mode
        weighted_gof_scores.append(score)
    
    return weighted_gof_scores

#Number of disjoint fROIs calculation
def number_disjoint_fROIs(zmap, threshold=3.09, n=1):
    num_disjoint_fROIs = []
    
    # Loop through each mode (4th dimension)
    for mode in range(zmap.shape[3]):  
        mode_data = zmap[:, :, :, mode]  # Extract the 3D volume for the current mode
        
        # Create a binary mask: True for voxels above threshold.
        activated = mode_data > threshold
        
        # Define the kernel size and build the cubic connectivity kernel.
        kernel_size = 2 * n + 1
        kernel = np.ones((kernel_size, kernel_size, kernel_size), dtype=int)
        kernel[n, n, n] = 0  # Exclude the center voxel from being counted.
        
        # Convolve the binary activation mask with the kernel to count activated neighbors.
        neighbor_counts = convolve(activated.astype(int), kernel, mode='constant', cval=0)
        
        # Identify isolated voxels: activated voxels with 0 activated neighbors.
        isolated_mask = activated & (neighbor_counts == 0)
        
        # Sum up the isolated voxels (each True counts as 1).
        num_disjoint_fROIs_each = np.sum(isolated_mask)
        num_disjoint_fROIs.append(num_disjoint_fROIs_each)
    
    return num_disjoint_fROIs

# Isoperimetric Quotient calculation
def iso_quotient(zmap, threshold=3.09):

    Q_values = []
    
    # Loop through each mode (4th dimension)
    for mode in range(zmap.shape[3]):  
        mode_data = zmap[:, :, :, mode]  # Extract the 3D volume for the current mode
        
        # Threshold the data to create a boolean activation mask
        mask = mode_data > threshold
        
        # Compute Volume (V) = total number of activated voxels
        V = np.sum(mask)
        if V == 0:
            Q_values.append(np.nan)
            continue
        
        # Compute Surface Area (S)
        # Use a 3x3x3 convolution to count the number of activated voxels in each neighborhood.
        kernel = np.ones((3, 3, 3), dtype=int)
        local_sum = convolve(mask.astype(int), kernel, mode='constant', cval=0)
        
        # A voxel is on the surface if it is activated and at least one neighbor is not activated (local_sum < 27).
        surface_mask = mask & (local_sum < 27)
        S = np.sum(surface_mask)
        if S == 0:
            Q_values.append(np.nan)
            continue
        
        # Compute Q using the formula:
        #    Q = ((4π)^(1.5) * V) / ((4/3)π * (S^(1.5)))
        numerator = (4 * math.pi) ** 1.5 * V
        denominator = ((4.0 / 3.0) * math.pi) * (S ** 1.5)
        Q = numerator / denominator
        
        Q_values.append(Q)
    
    return Q_values

#Fractal dimensions calculation
def compute_fractal_dimensions_4d(zmap, threshold=3.09, min_box_size=1, n_scales=10):
    fractal_dims = []
    num_modes = zmap.shape[3]
    sx, sy, sz = zmap.shape[:3]
    max_box_size = min(sx, sy, sz) // 2

    # Generate box sizes (scales) logarithmically spaced
    scales = np.unique(np.round(np.logspace(np.log10(min_box_size), np.log10(max_box_size), num=n_scales)).astype(int))
    
    for mode in range(num_modes):
        mode_data = zmap[..., mode]
        # Create a binary volume (activated if above threshold)
        binary_volume = (mode_data > threshold).astype(int)
        
        counts = []
        inv_scales = []
        # For each box size, count the number of boxes that contain at least one activated voxel.
        for box_size in scales:
            if box_size < 1:
                continue
            count = 0
            for i in range(0, sx, box_size):
                for j in range(0, sy, box_size):
                    for k in range(0, sz, box_size):
                        subvolume = binary_volume[i:min(i+box_size, sx),
                                                    j:min(j+box_size, sy),
                                                    k:min(k+box_size, sz)]
                        if np.any(subvolume):
                            count += 1
            counts.append(count)
            inv_scales.append(1.0/box_size)
        
        counts = np.array(counts)
        inv_scales = np.array(inv_scales)
        
        # Filter out scales with zero counts to avoid log(0)
        valid = counts > 0
        if np.sum(valid) < 2:
            # Not enough points to perform a linear regression
            fractal_dims.append(np.nan)
            continue

        log_inv_scales = np.log(inv_scales[valid])
        log_counts = np.log(counts[valid])
        
        # Fit a line to the log-log plot: log(N) vs. log(1/box_size)
        log_inv_scales = np.log(inv_scales)
        log_counts = np.log(counts)
        slope, intercept = np.polyfit(log_inv_scales, log_counts, 1)
        fractal_dims.append(slope)
    
    return fractal_dims

#Lacunarity calculation
def compute_lacunarity_4d(zmap, box_sizes, threshold=3.09, step=1):
    lacunarity_results = []
    sx, sy, sz, num_modes = zmap.shape

    for mode in range(num_modes):
        mode_data = zmap[..., mode]
        # Create a binary volume based on the threshold.
        binary_volume = (mode_data > threshold).astype(int)
        
        # Dictionary to store lacunarity values for each box size for this mode.
        lac_dict = {}
        
        for r in box_sizes:
            if r < 1:
                continue  # Skip invalid box sizes
            
            mass_count = {}  # Dictionary: mass (sum in a box) -> count
            box_count = 0    # Total number of boxes sampled at this scale
            
            # Slide a box of size r over the volume with given step.
            for i in range(0, sx - r + 1, step):
                for j in range(0, sy - r + 1, step):
                    for k in range(0, sz - r + 1, step):
                        subvol = binary_volume[i:i+r, j:j+r, k:k+r]
                        M = np.sum(subvol)  # mass in this subvolume
                        mass_count[M] = mass_count.get(M, 0) + 1
                        box_count += 1
            
            # If no boxes were sampled, set lacunarity to NaN.
            if box_count == 0 or len(mass_count) == 0:
                lac_dict[r] = np.nan
                continue
            
            # Compute the probability distribution Q(M, r) and then lacunarity.
            sum_M = 0.0
            sum_M2 = 0.0
            for M, nM in mass_count.items():
                Q_M_r = nM / box_count
                sum_M += M * Q_M_r
                sum_M2 += (M**2) * Q_M_r
            
            # Avoid divide-by-zero if sum_M is zero.
            if sum_M == 0:
                lac_dict[r] = np.nan
            else:
                L_r = sum_M2 / (sum_M ** 2)
                lac_dict[r] = L_r
        
        lacunarity_results.append(lac_dict)
    
    return lacunarity_results
