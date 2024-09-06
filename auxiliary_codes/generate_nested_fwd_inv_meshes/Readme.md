# MakeNestedMesh

**MakeNestedMesh** is a C++ program to generate the inversion and the nested forward modeling unstructured tetrahedral meshes, as well as the parameters of initial and/or reference models required by  **CRT3DMT** for inversion.

## 1 Prerequisites

### 1.1 Directory structure

The following files and sub-directories will be found:

**README.md**  Explains how to build and use this program.

**makefile** A file used by the `make` utility to automate the build process of this C++ program. It defines a set of rules to determine how to compile and link the program.

**MakeNestedMesh.cpp**  Entry of this program.

**em.h** Definition of namespace EM used in MakeNestedMesh.cpp .

**point.h** Definition of class  Point used in MakeNestedMesh.cpp.

**example/** Provides a example that is known to work.

### 1.2 Compiler

The Intel Compiler or GNU Compiler can be used.  The Intel Compiler can be accessed at link:  https://software.intel.com/en-us/intel-compilers. To use the Intel Compiler, `ENV{CXX}` in file `MakeNestedMesh/makefile` should be assigned to `icpc` and to used the GNU Compiler，that `ENV{CXX}` should be assigned to `g++` .

### 1.3 Third-party libraries

**TetGen**

Download the open source unstructured tetrahedral mesh generator TetGen (Version 1.4 is recommended). A copy of that library can be found in `CRT3DMT/contrib/tetgen1.4.3`.

### 1.3 Build

The building tool [GNU make](https://www.gnu.org/software/make/) should be installed before building the program. Steps to build the program:

(1) $ `cd` the `generate_nested_fwd_inv_meshes` directory

(2) $ `make clean`

(3) $ `make MakeNestedMesh`

After building the source codes, a executable program called **MakeNestedMesh** can be found in the generate_nested_fwd_inv_meshes` directory.

## 2 Usage

  The command to run the executable program  is              

```
./MakeNestedMesh [model_name] 
```

**The formats of input files and output files of this program are listed in the next section**.

## 3 File formats

In the following part, * denotes the model name.

### 3.1 Input files

 **(1) *.poly** 

*.poly file is the input model file required by the open source unstructured tetrahedral mesh generator **TetGen** (Version 1.4 is recommended) . For the detailed format of poly file, please refer to the User's Manual of TetGen (i.e., `CRT3DMT/contrib/tetgen1.4.3/manual1.4.pdf`).

**Please note **that for ***boundary marker***: use marker "2" for boundary of the computation domain, "1" for the air-earth-interface, and others for other boundary; for ***regional attribute***: use attribute "9999999" for the region of the air space, "6666666" for the non-active region in subsurface (if exists)  and other attributes for the other regions; For *node marker*: the markers  of observing sites must be different from the  markers of other nodes.

**(2)  elec.para**: 

elec.para is a file that describes the physical properties of different regions in the poly file. It has the following format:

```
    n_regions                           // the number of regions
    attribute sigma epsilon_r mu_r      // a list of attribute, conductivity (sigma)
                                        // relative dielectric constant and                                                       // relative magnetic permeability (mu_r)
                                        // for each 3-D region
```


Where, the words after double slash are the comments to the left content, which are not present in the file.

**(3) *.a.node**

*.a.node is a node file used for local refinement at the observing sites on the inversion mesh. The format of this file is the same as that of the node part in poly file.

**(4)*_inv.a.node**

*_inv.a.node is a node file used for local refinement at the observing sites on the nested forward modeling mesh. The format of this file is the same as that of the node part in poly file.

**For more details about the input files, please refer to the example files at the path `example\Inputs`**.

### 3.2 Output files

**(1) forward modeling and inversion mesh files** (.node, .ele, .face and .neigh)

For the detailed formats of the forward modeling mesh files (having keyword "fwd" in the file name, i.e. *_fwd.1.node), and the inversion mesh files (having keyword "inv" in the file name), please refer to the User's Manual of TetGen (i.e., `contrib/tetgen1.4.3/manual1.4.pdf`). 

**(2) *_inv.1.vol**:  

*_inv.1.vol is a volume constraint file used by refining the inversion mesh to generate the nested forward modeling mesh. Please refer to the User's Manual of TetGen (i.e., `contrib/tetgen1.4.3/manual1.4.pdf`) for the details about the format.

**(3) elec_regional_property.para**: 

The elec_regional_property.para file has the same format as that of the elec.para file. Which offers the parameters in the initial and reference model file .ini_ref required by **CRT3DMT** for inversion.

**(4) resistivity.vtk**

resistivity.vtk is a resistivity model file generated based on the regional resistivity information provided in the elec.para input file. Which can be displayed by using the visualization software  [ParaView](https://www.paraview.org/).

**For more details about the output files, please refer to the example files at the path `example\Outputs`**.