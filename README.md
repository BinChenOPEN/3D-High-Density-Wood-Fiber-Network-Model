# 3D High-density Wood Fiber Network Models
Author: Bin Chen

KTH Royal Institute of Technology, Sweden

E-mail 📧: binchen@kth.se or cbbuaa@outlook.com

## Project summary
This project aims to numerically generate 3D wood fiber network models with very high density (or low porosity). The 3D model can be used to simulate the physical properties of wood fiber network or paper materials. The advantage is the simplicity and high effeciency.


## 📔 Log
The software was confirmed to work on Windows 11, MATLAB 2020a.
***
## Parameters setting
**2025/09/27:**
1. Several parameters need to be set in ``Main.m``.
   - ``cell_wall_thick``, the thickness of the cell wall, in unit of voxel.
   - ``length_small``, the size of the middle model in x-y plane. The final volume will be slightly smaller than this value, here ``length_small-18`` voxels is final saved volume.
   -  ``fiber_align_mode``, indicate an isotropic or anisotropic in-plane fiber distributin by assigning value of 1 or 2.
   -  ``fiber_num``, How many fibers will be included, the more fibers meaning thicker material.
   -  ``fiber_width_mean``, the mean value of the fiber width.
   -  ``fiber_width_variation``, the variation of the fiber width in the range of ``[-fiber_width_variation,fiber_width_variation]``.
2. The obtained microstructure is saved in ``volum_compress_solid_center_new.mat``
3. The 3D microstructure can be viewed in the ``Volume Viewer`` app in Matlab
***
## Results
One microstructure with 50 fibers is shown in ``Figure 1``.


<p align="middle">
  <img src="Figure/Wood Fiber Network Isotropic 50.jpg" height="300" />
</p>

<p align="center"> Figure 1. 3D wood fiber network with 50 fibers. The fibers are isotropically distributed in plane.</p>


## Questions & Suggestions
Please contact Bin Chen (binchen@kth.se) if you wish to contribute code/algorithms to this project, or have question or suggestion. 

## Citation
Anyone who uses the code please cite: ***Simulated 3D microstructural geometries in wood fiber networks - towards high fiber content (to be submitted)***. If you need to redistribute the software, please keep the original author information.
