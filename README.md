# High-density Wood Fiber Network
Author: Bin Chen

KTH Royal Institute of Technology, Sweden

E-mail 📧: binchen@kth.se or cbbuaa@outlook.com

## Project summary
This project aims to generate wood fiber network with very high density (or low porosity). The 3D model can be used to simulate the physical properties of wood fiber network or paper materials. The model has following flexibility:

- Simple and effecient;
- Can control the fiber size and orientation distribution
- The model can reach a very low porosity (about 21%);
- Can simulated the molded shape, such as spherical and wedge shape;

## 📔 Log
**2025/09/27:**

The software was confirmed to work on Windows 11, MATLAB 2020a.


## Key information
### Parameters setting
1. Several parameters need to be set in ``Main.m``.
   - ``cell_wall_thick``, the thickness of the cell wall, in unit of voxel.
   - ``length_small``, the size of the middle model in x-y plane. The final volume will be slightly smaller than this value, here ``length_small-18`` voxels is final saved volume.
   -  ``fiber_align_mode``, indicate an isotropic or anisotropic in-plane fiber distributin by assigning value of 1 or 2.
   -  ``fiber_num``, How many fibers will be included, the more fibers meaning thicker material.
   -  ``fiber_width_mean``, the mean value of the fiber width.
   -  ``fiber_width_variation``, the variation of the fiber width in the range of ``[-fiber_width_variation,fiber_width_variation]``.
   - ``mold_type``, this paarameter should be 1 or 2 for spherical shape ro wedge shape
   - ``is_molded_WFN``, The value is selected between 0 and 1. 1 will generate the molded WFN with a complex shape, while 0 will skip that.
2. The obtained microstructure is saved in ``results/volum_compress_solid_center_new.mat``. While the compressed molded WFN microstructures are saved in ``results/volum_compress_curve_clean_center.mat``.
3. The 3D microstructure can be viewed in the ``Volume Viewer`` app in Matlab
3. If you want to rerun the project from a new set of parameters, please delete the corresponding result folder.

***
## Results
One microstructure with 50 fibers is shown in ``Figure 1``. The compressed high-density WFNs with spherical and wedge shapes are shown in ``Figure 2`` and ``Figure 3``, respectively.

<p align="middle">
  <img src="Figure/Wood Fiber Network Isotropic 50.jpg" width="300" />
</p>
<p align="center"> Figure 1. Compressed 3D wood fiber network with 50 fibers. The fibers are isotropically distributed in plane.</p>

<p align="middle">
  <img src="Figure/Wood Fiber Network Isotropic_spherical_shape_50.jpg" width="300" />
</p>
<p align="center"> Figure 2. Compressed 3D wood fiber network with 50 fibers and spherical shape. The fibers are isotropically distributed in plane.</p>

<p align="middle">
  <img src="Figure/Wood Fiber Network Isotropic_wedge_shape_50.jpg" width="300" />
</p>
<p align="center"> Figure 3. Compressed 3D wood fiber network with 50 fibers and wedge shape. The fibers are isotropically distributed in plane.</p>


## Questions & Suggestions
Please contact Bin Chen (binchen@kth.se), if you wish to contribute code/algorithms to this project, or have question or suggestion. 

## Citation
Anyone who uses the code please cite: ***Simulated 3D microstructural geometries in wood fiber networks - towards high fiber content (to be submitted)***. If you need to redistribute the software, please keep the original author information.
