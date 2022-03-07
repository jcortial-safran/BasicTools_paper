---
title: 'BasicTools: a numerical simulation toolbox'
tags:
  - Python
  - C++
  - mesh
  - fields
  - finite elements
  - post-treatment
authors:
  - name: Felipe Bordeu^[corresponding author]
    affiliation: 1
  - name: Fabien Casenave
    orcid: 0000-0002-8810-9128
    affiliation: 1
affiliations:
 - name: Safran Tech, Digital Sciences & Technologies Department, Rue des Jeunes Bois, Châteaufort, 78114 Magny-Les-Hameaux, France
   index: 1
date: 17 Février 2022
bibliography: paper.bib

---

# Summary

Numerical simulations of physical phenomena can be computed by many (commercially/free) software, but many similar operations are required when producing, preparing and consuming these simulations. Post-treatments of physical fields are a common need, as well as handling and modifying meshes and fields. BasicTools is a Python library designed especially to address these supporting tasks. It contains an efficient data model for meshes and field objects and input/output routines compatible with various formats. A finite element engine allows the assembling of abstract variational formulations, differential operators and the integration of fields on volume and surfaces.

BasicTools has been used in various projects in artificial intelligence and model order reduction for physical problems [@ROM-net; @mca26010017; @UQindustrialDesign; @datatargetVAE], topological optimization [@nardoni] and material sciences [@pymicro].

# Statement of need

When considering numerical simulation data from various sources, different formats may be handled. Each software can include its own post-treatment tool, which complexifies any comparison task or workflows mixing from different origins. Hence, common tasks and post-treatments may need to be reimplemented.

With BasicTools, we address these assessments by proposing a data model for meshes and associated physical fields, that is populated thanks to various reader and can be exported by various writers. Additionally, many treatment tools are implemented on these meshes and physical fields, including a finite element engine and field projection operators. We underline the fact that no new mesh or solution format is being specified: the goal is not to try and impose a new standard, but rather to handle existing data formats and propose common tasks and services on these data. Additionally BasicTools contain the tools to convert meshes to other "in memory" formats (Vtk [@VTK4], PyVista [@sullivan2019pyvista], MeshIO [@meshio]), this make possible to mix (and reuse) treatments already available in other frameworks.

# State of the field

In the computational fluid dynamics community, the CGNS [@cgns] standard has been proposed for analysis data. In the solid mechanics community, up to our knowledge, no such standard exists and a multiplicity of format exist. Concerning meshes, MeshIO has been proposed for conversion between various file format. For the manipulation of meshes Vtk is a good candidate but lacks some capabilities related to the solid mechanics community (e.g. integration point data). Finite element engines allowing to assemble abstract variational formulations are available in the codes FreeFem++[@freefempp] and FEniCS [@fenics].

# Overview

The architecture of the source code is described below.

```
BasicTools
└───Containers
│   │   ConstantRectilinearMesh.py
│   │   UnstructuredMesh.py
│   │   Filters.py
│   │   ...
│
└───FE
│   └───Fields
│   └───Integrators
│   └───Spaces
│   └───WeakForms
│   └───...
│   │
│   │   FETools.py
│   │   Integration.py
│   │   ...
│
└───IO
│   │   XdmfReader.py
│   │   XdmfWriter.py
│   │   ...
│
└───ImplicitGeometry
│   ...
```
</a>

The main objects and features of the library are

- meshes (in the folder Containers):
  ConstantRectilinearMesh.py and UnstructuredMesh.py define the data model for self-explicit mesh types. Unstructured meshes are efficient in the sense that elements are stored using only one array for each element type. Both mesh types can feature nodes and element tags. Many functions are available for creating, cleaning and modifying meshes. In particular, field projection operations enable to project fields defined on a mesh onto a set of points, using various methods and options, with respect to the location of the destination points being inside or outside the origin mesh (finite element interpolation, extrapolation, clamped evaluations, nearest neighbors, zero fill). Mesh morphing capacities are also included.
- filters (in the folder Containers):
  ElementFilter and NodeFilter enable to easily handle subparts of the meshes by selecting element and node set using threshold functions, tags, element types, element dimensionality and masks. This filter can be combined using boolean operations (union, complementary...) to construct advanced filters  on points and elements.
- a finite element engine (in the folder FE):
  A general weak formulation engine is available, enable to integrate fields over any part of the considered mesh. FETools.py contains particular functions for Lagrange P1 finite elements, including the computation of stiffness and mass matrices. The domain of integration is defined using the ElementFilters, this make the integration domain flexible. Depending on the parameter of the integration the result can be a matrix (i.e. tangent operator), vector (i.e. right hands side term), or a scalar (i.e. volume, energy)
- input/output functions (in the folder IO):
  Various readers (writer) for importing (exporting) meshes and solution fields from (to) the internal data model of BasicTools are available. Available formats include geo/geof (Z-set [@zset]), vtk, xdmf, samcef, abaqus. The bridge with MeshIO brings extra import/export capabilities by wrapping the MeshIO readers/writer with the same API as BasicTools.
- implicit geometry engine (in the folder ImplicitGeometry):
  This object are used to define geometries using only implicit geometries (level-set function). Basic shapes (spheres, plane, cylinders, cubes), operator (symmetry, translation, rotation) as well as binary operators (union, difference, intersection) can be used to construct complex shapes. Then these shapes can be used to select elements (using the ElementFilter), or be evaluated on a point cloud (e.g. the points of a mesh) to explicitly construct the iso-zero surface.

The large majority of functions are illustrated in the same file where they are defined, in CheckIntegrity functions. Additionally, some feature are illustrated below.

# Examples

We detail two examples illustrating some of the features listed in the previous section.
The complete examples can be found in the documentation of BasicTools.

## Pre/post deep learning

Some deep learning workflow applied to physics contexts require the projection of fields defined on an unstructured mesh onto a rectilinear grid, and inversely.

![Example of deep learning prepost.\label{fig:DeepLearningPrepost}](DeepLearningPrepost.png)

The mesh and solution associated to a previously computed flow field is read (see also \autoref{fig:DeepLearningPrepost}, top-left image):

```python
reader = XR.XdmfReader(filename = "PrePostDeepLearning_Input.xmf")
reader.Read()

grid = reader.xdmf.GetDomain(0).GetGrid(0)

# Read the mesh
uMesh = grid.GetSupport()
# convert all data to the correct binary (float64, int64) representation
# for the c++ part of BasicTools
uMesh.ConvertDataForNativeTreatment()

# Read the solution field 'U'
indexU = grid.GetPointFieldsNames().index("U")
U = grid.GetPointFields()[indexU][:,0:2]
```

A rectilinear mesh of size 48*48 excluding the left part of the mesh (x<0) is created:

```python
uMesh.ComputeBoundingBox()
Lx = uMesh.boundingMax[0] - 0.
Ly = uMesh.boundingMax[1] - uMesh.boundingMin[1]
rectMesh = ConstantRectilinearMesh(dim=2)
Nx = 48; Ny = 48
rectMesh.SetDimensions([Nx,Ny])
rectMesh.SetSpacing([Lx/(Nx-1), Ly/(Ny-1)])
rectMesh.SetOrigin([0.,uMesh.boundingMin[1]])
```

The projection operator from the unstructured mesh to the structured mesh is computed (see also \autoref{fig:DeepLearningPrepost}, top-right image):

```python
unstructuredRectMesh = CreateMeshFromConstantRectilinearMesh(rectMesh)
space, numberings, _offset, _nGauss =  PrepareFEComputation(uMesh)
inputFEField = FEField(name="U", mesh=uMesh, space=space, \
numbering=numberings[0])
method = "Interp/Clamp"
operator, status = GetFieldTransferOp(inputFEField, \
unstructuredRectMesh.nodes, method = method, verbose=True)

# Compute the projected field on the structured mesh
projectedU = operator.dot(U)
```

The projection operator from the structured mesh to the unstructured mesh (inverse projection) is computed (see also \autoref{fig:DeepLearningPrepost}, bottom-left image):

```python
space, numberings, _os, _nG = PrepareFEComputation(unstructuredRectMesh)
inputFEField = FEField(name="U", mesh=unstructuredRectMesh, \
space=space,numbering=numberings[0])
method = "Interp/Clamp"
operator, status = GetFieldTransferOp(inputFEField, uMesh.nodes, \
method = method, verbose=True)

# Compute the inverse-projected projected field on the unstructured mesh
inverseProjected_ProjectedU = operator.dot(projectedU)
```

The difference between the initial solution field and the result of the projection, then inverse projection of this field is illustrated on \autoref{fig:DeepLearningPrepost}, bottom-right image.



## Mechanical analysis: Thick plate with 2 inclusions

Here we present a study case of a thick plate with 2 inclusions. One softer and the other stiffer than the base material. We then compute the energy deformation on only one inclusion.

![Example of a mechanical thick plate with two inclusions analysis.\label{fig:TwoInclusions}](TwoInclusions_img1.png)

(We omitted the include statements for the sake of clarity)

```python
# Main class to drive a linear Finite Element analysis
problem = UnstructuredFeaSym()

# The mechanical problem (formulas for the weak form)
mecaPhys = MecaPhysics()
# Definition of the degree of the spaces [1 or 2]
mecaPhys.SetSpaceToLagrange(P=1)

# Add weak form terms to the tangent matrix
mecaPhys.AddBFormulation("Bulk",mecaPhys.GetBulkFormulation(1.0,0.3))
mecaPhys.AddBFormulation("Inclusion1",mecaPhys.GetBulkFormulation(5.0,0.3))
youngModulusInclusionII = 0.5
mecaPhys.AddBFormulation("Inclusion2",mecaPhys.GetBulkFormulation(0.5,0.3))

# Add weak form term to the rhs
mecaPhys.AddLFormulation( "Right", mecaPhys.GetForceFormulation([1,0,0],1))

# Push the physics into the FE problem
problem.physics.append(mecaPhys)

# The boundary conditions (block u,v,w on the left part of the mesh)
dirichlet = KRBlockVector()
dirichlet.AddArg("u").On('Left').Fix0().Fix1().Fix2()

problem.solver.constraints.AddConstraint(dirichlet)

# Read The mesh
mesh = ReadGmsh("TwoInclussions.msh")
mesh.ConvertDataForNativeTreatment()
print(mesh)

# Set the mesh into the problem
problem.SetMesh(mesh)

# We compute the numbering (mapping from the mesh to the linear system)
problem.ComputeDofNumbering()

# Assembly the linear system
with Timer("Assembly"):
    k,f = problem.GetLinearProblem()

# Compute the constraints to add to the system
problem.ComputeConstraintsEquations()

with Timer("Solve"):
    problem.Solve(k,f)

# Push the data from the solution vector to the FEFields
problem.PushSolutionVectorToUnkownFields()

# Recover a point representation of the displacement
# This is a no-op if the degree of the interpolation is 1 (P=1)
problem.mesh.nodeFields["sol"] = GetPointRepresentation(problem.unkownFields)

# Creation of a fake fields to export the rhs member
rhsFields = [FEField(mesh=mesh, space=None, \
numbering=problem.numberings[i]) for i in range(3)]
VectorToFEFieldsData(f,rhsFields)
problem.mesh.nodeFields["RHS"] = GetPointRepresentation(rhsFields)

print("Done solve")
print("Compute of the strain energy only on the second inclusion")
print("(integral in each element)")

symdep = GetField("u",3)
K = HookeIso(youngModulusInclusionII,0.3,dim=3)
# This is the projection space
symCellDataT = SWF.GetTestField("cellData",1)

print("Post process")
strainInVoigt = ToVoigtEpsilon(Strain(symdep))
EnerForm = strainInVoigt.T*K*strainInVoigt*symCellDataT

print("Post process Eval")
# Set the integration domain
ff = ElementFilter(mesh=problem.mesh, tag="Inclusion2")
p0Numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0)
energyDensityField = FEField(name="cellData", mesh=mesh, \
numbering=p0Numbering, space=LagrangeSpaceP0)

# This is the actual integration, because the target space is constant
# by element, we obtain for each element the strain energy
m, energyDensity = IntegrateGeneral(mesh=problem.mesh, wform=EnerForm,  \
constants={}, fields=problem.unkownFields, \
unkownFields = [energyDensityField], elementFilter=ff)

energyDensityField.data = energyDensity

problem.mesh.elemFields["Energy"] = energyDensity

import numpy as np
print("Strain energy on the second inclusion:", np.sum(energyDensity) )

writer = XdmfWriter.XdmfWriter('TwoInclussions_Output.xdmf')
writer.SetHdf5(False)
writer.Open()
writer.Write(mesh,PointFields=list(mesh.nodeFields.values()), \
                  PointFieldsNames=list(mesh.nodeFields.keys()), \
                  CellFields=list(mesh.elemFields.values()), \
                  CellFieldsNames=list(mesh.elemFields.keys()))
writer.Close()
```

# Acknowledgements

The authors wish to thank Julien Cortial (Safran) for his contributions to the library.

# References