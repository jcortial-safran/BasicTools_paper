from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh
from BasicTools.Containers.UnstructuredMeshModificationTools import Morphing, ComputeSkin
from BasicTools.Containers.Filters import NodeFilter, ElementFilter 
from BasicTools.IO import XdmfWriter as XW
import numpy as np


############################
# CONSTRUCT AN INITAL MESH #
############################

# Construct a rectilinear mesh of size 8*12*16 
structuredMesh = ConstantRectilinearMesh(dim=3)
structuredMesh.SetDimensions([8,12,16])
structuredMesh.SetOrigin([0.,0.,0.])
structuredMesh.SetSpacing([0.1,0.1,0.1])


# Construct a rectilinear mesh of size 8*12*16 
unstructuredMesh = CreateMeshFromConstantRectilinearMesh(structuredMesh,ofTetras=True)


# Export the unstructured mesh in xdmf format
writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='init.xdmf')
writer.Write(unstructuredMesh, PointFields = [], PointFieldsNames = [])
writer.Close()


##################
# MORPH THE MESH #
##################

# Construct an element tag corresponding to the skin of the mesh 
ComputeSkin(unstructuredMesh, inplace=True)


# Construct an nodeFilter filtering nodes of the skin such that z>1
nodeFilter = NodeFilter(unstructuredMesh)
nodeFilter.AddETag("Skin")
nodeFilter.AddZone(lambda p: 1-p[:,2])

# Add a nodeTag tagging the elements of nodeFilter
idsToTreat = nodeFilter.GetIdsToTreat()
nodalTag = unstructuredMesh.GetNodalTag("ImposedDisp")
nodalTag.AddToTag(idsToTreat)


# compute a displacement to apply to the selected nodes
posIdsToTreat = unstructuredMesh.GetPosOfNodes()[idsToTreat]

dispIdsToTreat = np.zeros(posIdsToTreat.shape)
dispIdsToTreat[:,0] = (posIdsToTreat[:,2] - 1.) * (posIdsToTreat[:,0] - 1.)
dispIdsToTreat[:,1] = (posIdsToTreat[:,2] - 1.) * (posIdsToTreat[:,1] - 1.)


# morph the mesh by updated the position of the nodes
morphedNodes = Morphing(unstructuredMesh, dispIdsToTreat, idsToTreat)
unstructuredMesh.SetNodes(morphedNodes)


###############################################
# PROJECT A SCALAR FIELD ON THE MORPHDED MESH #
###############################################

# Construct a rectilinear mesh of size 10*10*10 for the bounding box of the morphed mesh
unstructuredMesh.ComputeBoundingBox()

structuredMeshBB = ConstantRectilinearMesh(dim=3)
n = 10
spacing = [(unstructuredMesh.boundingMax[i]-unstructuredMesh.boundingMin[i])/(n-1) for i in range(3)]
structuredMeshBB.SetDimensions([n,n,n])
structuredMeshBB.SetOrigin(unstructuredMesh.boundingMin)
structuredMeshBB.SetSpacing(spacing)


# Construct a field supported on this rectilinear mesh
testField = np.empty(n**3)
structuredMeshBBNodes = np.empty((n**3,3))
count = 0
for i in range(n):
    x = i*spacing[0]
    for j in range(n):
        y = j*spacing[1]
        for k in range(n):
            z = k*spacing[2]
            testField[count] = np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)
            structuredMeshBBNodes[count,:] = [x,y,z]
            count += 1


# Export the rectilinear mesh and constructed field in xdmf format
writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='boundingBox.xdmf')
writer.Write(structuredMeshBB, PointFields = [testField], PointFieldsNames = ["testField"])
writer.Close()


# Compute the projection operator from the rectilinear mesh to the morphed unstructured mesh
from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh
from BasicTools.FE.FETools import PrepareFEComputation
from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
from BasicTools.FE.Fields.FEField import FEField

unstructuredMeshBB = CreateMeshFromConstantRectilinearMesh(structuredMeshBB)
space, numberings, _offset, _NGauss = PrepareFEComputation(unstructuredMeshBB, numberOfComponents = 1)
inputFEField = FEField(name="testField",mesh=unstructuredMeshBB,space=space,numbering=numberings[0])
operator, status = GetFieldTransferOp(inputFEField, unstructuredMesh.nodes, method = "Interp/Clamp", verbose=True)


# Compute the projected field on the morphed unstructured mesh
projectedTestField = operator.dot(testField)


# Export the morphed unstructured mesh and projected field in xdmf format
writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='morphed.xdmf')
writer.Write(unstructuredMesh, PointFields = [projectedTestField], PointFieldsNames = ["projectedTestField"])
writer.Close()


##############################################################
# INTEGRATE THE SCALAR FIELD ON THE SKIN OF THE MORPHED MESH #
##############################################################

# Construct an nodeFilter filtering elements of the skin
elFilter = ElementFilter(unstructuredMesh)
elFilter.AddTag("Skin")

from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.Helpers.Timer import Timer

print("--")
string = "Integrate the scalar field on the skin of the morphed mesh by applying numerical quadrature\n"
TFormat.II(2)
string += TFormat.GetIndent() + " from Lagrange P1 finite element integration using two different methods"
print(string)
TFormat.Reset()

#### Method 1
print(TFormat.Center(TFormat.InRed("Method 1:")+TFormat.InBlue("'by hand'")))

timer = Timer("Duration of method 1").Start()

from BasicTools.FE.FETools import ComputePhiAtIntegPointFromElFilter, ComputeL2ScalarProducMatrix

integrationWeights, phiAtIntegPoint = ComputePhiAtIntegPointFromElFilter(unstructuredMesh, elFilter)
integrationResult1 = np.dot(integrationWeights, phiAtIntegPoint.dot(projectedTestField))

timer.Stop()

print("integration result =", integrationResult1)

#### Method 2
print(TFormat.Center(TFormat.InRed("Method 2:")+TFormat.InBlue("using the weak form engine")))

timer = Timer("Duration of method 2").Start()

from BasicTools.FE.Integration import IntegrateGeneral
from BasicTools.FE.SymWeakForm import GetField, GetTestField
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1, ConstantSpaceGlobal
from BasicTools.FE.DofNumbering import ComputeDofNumbering

F = GetField("F",1)
Tt = GetTestField("T",1)

wf = F.T*Tt

numbering = ComputeDofNumbering(unstructuredMesh,LagrangeSpaceP1)
field = FEField("F",mesh=unstructuredMesh,space=space,numbering=numbering, data = projectedTestField)

numbering = ComputeDofNumbering(unstructuredMesh,ConstantSpaceGlobal)
unkownField = FEField("T",mesh=unstructuredMesh,space=ConstantSpaceGlobal,numbering=numbering)

K, F = IntegrateGeneral(mesh=unstructuredMesh,
                    wform=wf,
                    constants={},
                    fields=[field],
                    unkownFields=[unkownField],
                    integrationRuleName="LagrangeP1",
                    elementFilter=elFilter)

integrationResult2 = F[0]
timer.Stop()

print("integration result =", integrationResult2)
############

print(Timer.PrintTimes())



