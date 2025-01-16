# CRT3DMT

[![DOI](https://zenodo.org/badge/853193038.svg)](https://doi.org/10.5281/zenodo.13711326)

**CRT3DMT** is a three-dimensional (3-D) magnetotelluric (MT) forward modeling and inversion program using unstructured tetrahedral meshes. In terms of unstructured finite element solver, it can deal with arbitrarily complicated 3-D models, including rugged topographies. To enhance the reliability of inversion results and reduce computational costs, the nested forward modeling and inversion tetrahedral meshes along with an adaptive refinement strategy for the inversion mesh can be employed. You can use it for

* Producing 3-D MT synthetic impedance and tipper data with and without Gaussian noise.
* 3-D MT inversion using same forward modeling and inversion meshes.
* 3-D MT inversion using fixed nested forward modeling and inversion meshes.
* **3-D MT adaptive inversion using nested forward modeling and inversion meshes**.

## 1 Prerequisites

### 1.1 Directory structure

The following files and sub-directories will be found:

**LICENSE**  Please get familiar with the License before using this program!

**makefile** A file used by the `make` utility to automate the build process of this C++ program. It defines a set of rules to determine how to compile and link the program.

**README.md**  Explains how to build and use this program.

**CRT3DMT.cpp**  Entry of this program.

**src/**  Others c++ source files of this program.

**contrib/**  Open source packages that the program depends on directly or indirectly.

**auxiliary_codes/** Include the program to generate inversion mesh and the nested forward modeling mesh

**examples/** Examples for 3-D MT forward modeling and inversion.

### 1.2 Compiler

Install the free-version Intel Compiler because this program uses its PARDISO solver in Math Kernel Library.  The Intel Compiler can be accessed at link:  https://software.intel.com/en-us/intel-compilers. To use the Intel Compiler and PARDISO solver, the `ENV{CXX}` and  `ENV{PARDISO_LIB}` in file `CRT3DMT/makefile` should be assigned to `icpc` and `-mkl`, respectively. To enable OpenMP parallelization over observing frequencies, the flag `-qopenmp` should be added into `ENV{CXXFLAGS}` in file  `CRT3DMT/makefile` .

### 1.3 Third-party libraries

 **(1) Eigen**

Download the open source library Eigen used for basic linear algegra operations at link: http://eigen.tuxfamily.org/ and put it into folder `CRT3DMT/contrib`.  A copy of that library can be found in `CRT3DMT/contrib/eigen3.3.9`. To use Eigen library, the `ENV{EIGEN_INCLUDE}` in file `CRT3DMT/makefile` should be assigned to its path.

**(2) TetGen**

Download the open source unstructured tetrahedral mesh generator TetGen (Version 1.4 is recommended), which is used for adaptive inversion mesh refinement. A copy of that library can be found in `CRT3DMT/contrib/tetgen1.4.3`.

### 1.4 Build

The building tool [GNU make](https://www.gnu.org/software/make/) should be installed before building the program. Steps to build the program:

(1) $ `cd` the `CRT3DMT` directory

(2) $ `make clean`

(3) $ `make CRT3DMT`

After building the source codes, a executable program called **CRT3DMT** can be found in the `CRT3DMT` directory.

## 2 Usage

  The unified command to run the executable program for both forward modeling and inversion is

```
./CRT3DMT [startup_file_name]
```

The startup file contains two columns: the left columns are comments to the content on the right column. The comments at the left column cannot contain any space, as it is read as a string. 

**For forward modeling option, the startup file has the following format**:

```
Project_name:                              adinvem
Data_method:                               MT3D
Function:                                  Fwd_only
Forward_modeling_parameter_file_name:      crt3dmt_input_fwd.para
Rel_impedance_error(usually_in_0-0.1):     0.05
Abs_impedance_error(depends_on_abs(d)):    0.00
Rel_tipper_error(commonly_0.0):            0.00
Abs_tipper_error(commonly_0.01-0.03):      0.005
```

In which, the string "adinvem" denotes the project name, mandatory strings "Fwd_only" and "MT3D "  respectively denotes performing forward modeling and generating  3-D MT synthetic data, "crt3dmt_input_fwd.para" is the name of the forward modeling parameter file, the numbers  at the last four lines on right columns are used to contaminate the impedance and tipper modeling data.

**For inversion option, the startup file has the following format**:

```
Project_name:                              adinvem
Data_method:                               MT3D
Function:                                  L-BFGS
Data_file_name:                            crt3dmt_input.data
Initial_priori_model_file_name:            crt3dmt_input.ini_ref
Inversion_parameter_file_name:             crt3dmt_input.inv 
Forward_modeling_control_file_name:        crt3dmt_input.fwd 
```

In which, the string "L-BFGS" denotes performing L-BFGS inversion, "crt3dmt_input.data" is the name of the data file, “crt3dmt_input.ini_ref” is the name of the file contains initial and priori models, "crt3dmt_input.inv " is the name of inversion parameter file, and "crt3dmt_input.fwd" is the name of parameter contol file for forward modeling used in each L-BFGS iteration.

**The formats of other input files and output files for both forward modeling and inversion are listed in the next section**.

## 3 File formats

### 3.1 Input files

#### 3.1.1 forward modeling

 **(1) forward modeling parameter file** (*.para)

The forward modeling parameter file has the following format:

```
    n_frequencies                  // total number of observing frequencies
    f_1  f_2 ... f_n               // specific observing frequencies (Hz)
    mesh_1 mesh_2 ... mesh_n       // mesh name = mesh_i.counter.* for frequency f_i
    counter                        // * are node, ele, face and neigh.
    theta                          // incidence angle of electromagnetic wave,
                                   // theta = 0 denotes vertical incidence
    n_layers                       // total number of layers for 1-D background model
    sigma mu_r epsilon_r  d        // a list of conductivity (sigma), 
                                   // relative magnetic permeability (mu_r),
                                   // relative dielectric constant (epsilon_r),
      			                       // and depth (d) from air-earth interface (Z=0) 
                                   // to the bottom surface of each layer 
    n_regions                      // the number of regions
    attribute sigma mu_r epsilon_r // a list of attribute, conductivity (sigma)
                                   // relative mu (mu_r) and relative epsilon 
                                      (epsilon_r) for each 3-D region
    site_file_name                 // the name of site file
```

Where, the words after double slash are the comments to the left content, which are not presented in the file.

**(2) site file**

The site file has the following format:

```
    int N                // total number of observing sites  
    int marker           // site marker                      
    u_1                  // filds along u direction at site 1
    ...                  // filds along u direction at other N-2 sites
    u_N                  // filds along u direction at site N
    v_1                  // filds along v direction at site 1
    ...                  // filds along v direction at other N-2 sites
    v_N                  // filds along v direction at site N
    mu_r_1               // relative mu at site 1 for calculating impedance
    ...                  // relative mu at other (N-2) sites  
    mu_r_N               // relative mu at site N
```

In which, the words after double slash are the comments to the left content.

**(3) mesh files** (*.ele, *.node, *.neigh and *.face)

For the formats of TetGen mesh files: *.ele, *.node, *.neigh and .face , please refer to the User's Manual of TetGen (i.e., `contrib/tetgen1.4.3/manual1.4.pdf`). 

**For more details about the input files of forward modeling, please refer to the example at the path `examples\forward_modeling\Inputs`**.

#### 3.1.2 Inversion

**(1) data file** (*.data)

To maintain consistency with the data format used by ModEM, the data format adopted by this program has only the following differences compared with ModEM's data format: 

* Only the following four data type keywords are used:

  | Keyword                            | Allowed comonents          | Real or Complex |
  | :--------------------------------- | -------------------------- | --------------- |
  | Full_impedance                     | Zxx, Zxy, Zyx, Zyy         | Complex         |
  | Off_Diagonal_Impedance             | Zxy, Zyx                   | Complex         |
  | Full_Impedance_plus_Tipper         | Zxx, Zxy, Zyx, Zyy, Tx, Ty | Complex         |
  | Off_Diagonal_Impedance_plus_Tipper | Zxy, Zyx, Tx, Ty           | Complex         |

* For data type keyword "Full_Impedance_plus_Tipper" or "Off_Diagonal_Impedance_plus_Tipper", the tipper data must be placed right after the impedance data without head lines.

* Both impedance and tipper data must be stored in the following order: outer cycle is for site, middle cycle is for frequency (the order of frequencies for each station should be consistent), and  inner cycle is for data component (impedance component are stored in the order Zxx (if needed), Zxy, Zyx, Zyy (if needed)；tipper component  (if needed) are stored in the order (Tx, Ty)).

* The order of sites in  the data file must be consistent with the order of that in the mesh node file (*.node), as the program tracks sites based on that order in the *.node file.

Note: the coordinates of the sites provided in the data file are only used to check if their order is consistent with the order of that in the mesh node file and do not participate in the calculation. The program identifies the site location information by recognizing the marker information in the mesh node file.  

For more details about the data file, please refer to section  “4.1  data files” in document at the path`contrib/ModEM_UserGuide2019.pdf` . 

**(2) initial and reference model file** (*.ini_ref)

The initial and reference model file has the following format:

```
#starting_model,counter,n_regions(!Don'T_TYPE_SPACE_IN_THIS_LINE)       
half_space_two_layer_sites_inv
1
263673

#sigma,epsilon_r,mu_r(for_initial_model)(!Don'T_TYPE_SPACE_IN_THIS_LINE)
1	0.01	0	1	
2	0.01	0	1	
3	0.01	0	1	
...

#sigma,epsilon_r,mu_r(for_reference/prior_model)(!Don'T_TYPE_SPACE_IN_THIS_LINE)
default
```

As shown above, the content of this file consists of three parts. Each part begins with a comment line, which will be read as a string. 

**For the first part**: The second line and third line (counter) constitute the inversion mesh name. In this example, the inversion mesh name is half_space_two_layer_sites_inv.1.*，in which * are node, ele, face and neigh. The fourth line is the number of the regions in the inversion mesh, including active and fixed regions. Each regions has a unique regional attribute.

**For the second part**: This part describes the initial model. Each line except the first line shows the regional attribute, conductivity (sigma), relative dielectric constant (epsilon_r) and relative magnetic permeability (mu_r) of each 3-D region. 

**For the third part**: This part describes the reference model. The reference model can be described in the same way as that used for the initial model, and can also be described by using a keyword "default" which means using the initial model as the reference model.

**(3) inversion parameter file** (*.inv)

The inversion parameter file has the following format:

```
Adjustment_approach_of_inv_mesh:           1              
Max_No_of_inversion_unknowns:              10000000
Max_adjust_times_of_inv_mesh:              6
Max_inv_iter_times_per_mesh:               6
Frac_of_grad_m_to_refine(0-1.0):	        0.2 
Frac_of_others_to_refine(0-1.0):           0.2

Cond_lower_bound:                          1.0E-5
Cond_upper_bound:                          1.0E+1
Constant_for_cond_transformation(>0):      1.0
Lambda0:                                   100
Tol_reduc_frac_of_rms_to_cool_lambda:      0.2
Cool_factor:                               3
Tol_times_of_continuously_cool_lambda:     10
Tol_lambda:                                1.0E-14
Max_iteration_No:                          100
Tol_RMS:                                   1.0E-8
```

The inversion parameter file contains two columns: the left columns are comments to the parameter on the right column. The comments at the left column cannot contain any space, as it is read as a string. 

**1st line**:  Refining inversion mesh or not. "1" denotes refining inversion mesh adaptively (adaptive inversion), "0" denotes using fixed inversion mesh (traditional inversion). When choosing "0", the parameters for inversion mesh refinement at the second to sixth lines will be ignored.

**2nd ~  6th lines**: Parameters for inversion mesh refinement.  *Max_No_of_inversion_unknowns*: maximum allowable number of inversion unknowns; *Max_adjust_times_of_inv_mesh*: maximum allowable refinement times of inversion mesh ; *Max_inv_iter_times_per_mesh*: maximum allowable iteration times for each inversion mesh; *Frac_of_grad_m_to_refine*: refinement ration of inversion elements for refinement indicator of model gradient; *Frac_of_others_to_refine*: refinement ration of inversion elements for refinement indicator of data misfit gradient .

**7th ~ 9th lines**: Parameters for conductivity. *Cond_lower_bound* and *Cond_upper_bound* denote the lower and upper bound of conductivity, respectively.  *Constant_for_cond_transformation*  denotes the type of logarithmic transformation for conductivity, '1' denotes using natural logarithm, '2.30' denotes using common logarithm.

**10th ~ 12th lines**: Parameters for regularization factors. *Lambda0* is the initial regularization factor.*Tol_reduc_frac_of_rms_to_cool_lambda* is cooling threshold. When the decrease ratio of RMS falls below this cooling threshold，the regularization factor is devided by *Cool_factor*.

**13th ~ 16th lines**: Parameters for stopping criterion of inversion. *Tol_times_of_continuously_cool_lambda*: the maximum allowable consecutively decreased times of regularization factor; *Tol_lambda*: the minimum  allowable value of regularization factor; *Max_iteration_No*: maximum allowable times of inversion iterations; *TOl_RMS*: target value of RMS.

**(4) forward modeling control file** (*.fwd)

The forward modeling control file has the following format:

```
#starting_model,starting_counter,theta(!Don'T_TYPE_SPACE_IN_THIS_LINE)
half_space_two_layer_sites_fwd
1
0.0

#n_layers,sigma,epsilon_r,mu_r,d(used_for_boundary_condition)(!Don'T_TYPE_SPACE_IN_THIS_LINE)
2
1E-16  0.0    1.0    0.0
1E-02  0.0    1.0    1E50

#file_name_for_sites(!NO_SPACE)
input_file.site
```

As shown above, the content of this file consists of three parts. Each part begins with a comment line, which will be read as a string. 

**For the first part**: The second line and third line (counter) constitute the name of the forward modeling mesh . In this example, the name of the forward modeling mesh is half_space_two_layer_sites_fwd.1.*，in which * are node, ele, face and neigh. The fourth line is  the incidence angle of electromagnetic wave, theta = 0 denotes vertical incidence.

**For the second part**: This part describes 1-D background model used for loading boundary condition. The first line is the number of layers. Each successive  line shows the conductivity, relative dielectric constant, relative magnetic permeability, and depth from air-earth interface (z=0) to the bottom surface of each layer.

**For the third part**: The name of the site file.

**(5) site file** 

Using the same format as that of the site file for forward modeling, please refer to subsection 3.1.1 for more details.

**(6) forward modeling and inversion meshes**  (*.ele, *.node, *.neigh and *.face)

Please refer to the User's Manual of TetGen (i.e., `contrib/tetgen1.4.3/manual1.4.pdf`) for details of the formats of mesh files: *.ele, *.node, *.neigh and *.face.

**Note**: As each inversion element needs a unique attribute, the inversion mesh files should be generated by using the **MakeNestedMesh** program at the path `auxiliary_codes/generate_fwd_inv_meshes`. As the parameter transfer  between forward modeling mesh and inversion mesh is based on the attribute, the forward modeling mesh files can be the same as the inversion mesh files or generated by using **MakeNestedMesh** program as well.

**For more details about the input files of inversion, please refer to "Inputs" folder under the path of `examples\inversion`**.

### 3.2 Output files

#### 3.2.1 forward modeling

**(1) MT3D_Synthetic.data**: file storing synthetic data generated by using unstructured FEM

The format of synthetic data file is the same as that of data file for inversion. Please refer to subsection “3.1.2 Inversion“ for the details.

**(2) input_parameter_list.log**: log file for input parameter, used for double checking.

**(3) screen_output.log**：log file for screen output.

**For more details about the input files of forward modeling, please refer to the example at the path "\examples\forward_modeling\Outputs"**.

#### 3.2.2 Inversion

**(1) .RMS file** : storing the RMS values for each L-BFGS iteration.

**(2) .lambda file**: storing the value of regularization factor for each L-BFGS iteration.

**(3) .roughness file**: storing the roughness of the updated model for each L-BFGS iteration.

**(4) .search_times file**: storing the search times of inexact linear search for each in L-BFGS iteration.

**(5) .alpha file**: storing the step length of each L-BFGS iteration.

**(6) .n.rho.vtk file**: storing the resistivity model obtained at the n-th L-BFGS iteration.

**(7) .n.dat file**: storing the predicted data at the n-th L-BFGS iteration.

**(8) .n.res file**: storing the absolute deviation between the input data and the predicted data at the the n-th L-BFGS iteration.

**(9) input_parameter_list.log**: log file for input parameter, used for double checking.

**(10) screen_output.log**: log file for screen output.

**Note**: For adaptive inversion, the number of inversion mesh appears in the name of the files (1)-(8). i.e., ***.m.RMS** and ***.m.n.rho.vtk**, which denote storing the RMS values for each L-BFGS iteration on m-th inversion mesh and the resistivity model obtained at the n-th L-BFGS iteration on the m-th inversion mesh, respectively. Moreover, the refined inversion mesh files and the corresponding forward modeling mesh files are generated for adaptive inversion.

**For more details about the output files of inversion, please refer to "Outputs" folder under the path of `examples\inversion`**. 
