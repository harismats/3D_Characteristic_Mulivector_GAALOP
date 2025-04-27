# Code for Paper "Symbolically Optimized Characteristic Multivector Rotor Estimation for 3D Point Cloud Registration"

## Overview
This repository benchmarks Geometric Algebra–based rotation estimation against standard SVD within both Absolute Orientation (AO) setting and an Iterative Closest Point (ICP) framework. We:

1. Apply three Characteristic Multivector (CM) methods—  
   - **Original CM** (Clifford Toolbox)  
   - **CM GAALOP** (symbolically optimized)  
   - **CM GAALOP CSE** (with common‐subexpression elimination)  
   —to a sample point cloud (`apple.ply, horse.ply`) and compare runtime & accuracy against ground truth.

2. Embed those four methods (plus **SVD**) in ICP loops over synthetic clouds of increasing size to compare performance and error scaling.

## Requirements

- **MATLAB** R2020a or later  
- **Clifford Multivector Toolbox** (add with `addpath`)  
- **Computer Vision Toolbox** (for `pcread` & `pointCloud`)  
- Sample point cloud: `apple.ply`

## Setup

1. In MATLAB, add the toolbox and your functions folder to the path:
   ```matlab
   addpath(genpath('path/to/CliffordToolbox'))
   addpath(pwd)
   ```
2. Ensure `apple.ply, horse.ply` are in your current folder (or update its path in the scripts).

## Running the Benchmarks

### 1. Absolute Orientation (`main_AO.m`)

- **Purpose**: Compare three CM variants on 200 random rigid transforms.  
- **Outputs**:  
  - **Runtime vs. Test Case** (linear & log scale)  
  - **Rotation Error vs. Test Case**  
- **Usage**:
  ```matlab
  >> main_AO
  ```
- **Adjustable parameters** (at top of `main_AO.m`):  
  - `no_test_cases` – number of random transforms  
  - `outputFolder` – where figures are saved  
  - Point cloud filename

### 2. ICP Comparison (`main_ICP.m`)

- **Purpose**: Measure runtime & Frobenius‐norm rotation error for SVD, CM, CM GAALOP, CM GAALOP CSE across cloud sizes.  
- **Settings**:
  - Cloud sizes: 1000 : 500 : 8000 points  
  - 5 trials per size  
  - 100 ICP iterations each  
- **Outputs**:  
  - **Avg. Execution Time vs. Number of Points**  
  - **Avg. Rotation Error vs. Number of Points**  
  - **Percent‐Difference Plots** (relative to SVD)
- **Usage**:
  ```matlab
  >> main_ICP
  ```
- **Adjustable parameters** (in `main_ICP.m`):  
  - `pointCloudSizes`, `noTestCases`, `noIteration`  
  - `outputFolder`

## Core Concepts

Although the Characteristic Multivector (CM) method in Geometric Algebra yields high accuracy, its direct MATLAB implementation is relatively slow. By exporting its symbolic rotor formulas through GAALOP—and further applying CSE we generate two optimized variants that:

- Retain CM’s high accuracy  
- Reduce computation time
- Perform on par with SVD inside an ICP loop  

This demonstrates GAALOP‐optimized GA as a practical alternative for high‐accuracy 3D registration.

## Contact

**Haris Matsantonis**  
Email: cm2128@cam.ac.uk  

---
