from BasicTools.IO import XdmfReader as XR
import numpy as np

reader = XR.XdmfReader(filename = "jet_flame.sol000829.xmf")
reader.Read();
dom = reader.xdmf.GetDomain(0);

grid = dom.GetGrid(0)
unstructuredMesh = grid.GetSupport()
Uexact = grid.GetPointFields()[0][:,0:2]



#Read projected mesh and fields
from BasicTools.IO import XdmfReader as XR
import numpy as np

reader = XR.XdmfReader(filename = 'TestProj.xmf')
reader.Read();
dom = reader.xdmf.GetDomain(0);


grid = dom.GetGrid(0)

PointFieldsNames = grid.GetPointFieldsNames()
U = np.array(grid.GetPointFields())



Nx = 48
Ny = 48

originMesh = [0., np.min(unstructuredMesh.nodes[:,1])]
coordMax = [np.max(unstructuredMesh.nodes[:,0]), np.max(unstructuredMesh.nodes[:,1])]

Lx = coordMax[0] - originMesh[0]  
Ly = coordMax[1] - originMesh[1]  

#Create Projected mesh  
from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh    
rectMesh = ConstantRectilinearMesh(dim=2)
rectMesh.SetDimensions([Nx,Ny])
rectMesh.SetSpacing([Lx/(Nx-1), Ly/(Ny-1)])
rectMesh.SetOrigin([originMesh[0],originMesh[1]])


import BasicTools.Containers.UnstructuredMeshCreationTools as UMCT
cartesianMeshUnstructured = UMCT.CreateMeshFromConstantRectilinearMesh(rectMesh)


from BasicTools.FE.FETools import PrepareFEComputation
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
space, numberings, _offset, _NGauss = PrepareFEComputation(cartesianMeshUnstructured)
inputFEField = FEField(name="U",mesh=cartesianMeshUnstructured,space=space,numbering=numberings[0])
methods = ["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"]
method = methods[2]



operator, status = GetFieldTransferOp(inputFEField, unstructuredMesh.nodes, method = method, verbose=True)

print("operator.shape =", operator.shape)
print("U.shape =", U.shape)

InvProjFields = operator.dot(U.T).T

print("InvProjFields.shape =", InvProjFields.shape)



from BasicTools.IO import XdmfWriter as XW
#write solution in xdmf
writer = XW.XdmfWriter('TestInverseProj.xmf')
writer.SetTemporal(True)
writer.SetHdf5(False)
writer.Open()
writer.Write(unstructuredMesh,PointFields=[InvProjFields[i] for i in range(len(PointFieldsNames))], PointFieldsNames=PointFieldsNames)
writer.Close()






"""
#construct mesh for combustor zone
from BasicTools.Containers.Filters import ElementFilter
ff = ElementFilter(unstructuredMesh, zone = lambda p: (p[:,0]),zones=[])
plaqueElements = {}
for name,data,ids in ff:
    plaqueElements[name] = ids


from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes

#Read unstructuedMesh for current solution
reader = VR.VtuReader(nameFileInitSol)
meshin = reader.Read()
unstructuredMesh = ExtractElementsUsingDict(meshin,plaqueElements)
CleanLonelyNodes(unstructuredMesh)
"""





#calcul erreur
print("calcul erreur")
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.FE import FETools as F
from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByElementFilter
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes
from BasicTools.Containers.UnstructuredMeshFieldOperations import CopyFieldsFromOriginalMeshToTargetMesh


unstructuredMesh.nodeFields['Uexact'] = Uexact
unstructuredMesh.nodeFields['UInvProj'] = InvProjFields.T
unstructuredMesh.nodeFields['delta'] = Uexact-InvProjFields.T
ff = ElementFilter(unstructuredMesh, zone = lambda p: (-p[:,0]))



unstructuredMeshClipped = ExtractElementsByElementFilter(unstructuredMesh,ff)
CleanLonelyNodes(unstructuredMeshClipped)
CopyFieldsFromOriginalMeshToTargetMesh(unstructuredMesh,unstructuredMeshClipped)

#ff = ElementFilter(mesh)
#M = F.ComputeL2ScalarProducMatrix(unstructuredMesh, 2)#, elementFilter = ff)
integrationWeights, phiAtIntegPoint = F.ComputePhiAtIntegPoint(unstructuredMeshClipped)
print("integrationWeights.shape =", integrationWeights.shape)
print("phiAtIntegPoint.shape =", phiAtIntegPoint.shape)
print("unstructuredMeshClipped.nodeFields['delta'].shape =", unstructuredMeshClipped.nodeFields['delta'].shape)


vectDeltaAtIntegPoints = np.empty((2,phiAtIntegPoint.shape[0]))
vectUexactAtIntegPoints = np.empty((2,phiAtIntegPoint.shape[0]))
vectUInvProjAtIntegPoints = np.empty((2,phiAtIntegPoint.shape[0]))
for i in range(2):
    vectDeltaAtIntegPoints[i] = phiAtIntegPoint.dot(unstructuredMeshClipped.nodeFields['delta'][:,i])
    vectUexactAtIntegPoints[i] = phiAtIntegPoint.dot(unstructuredMeshClipped.nodeFields['Uexact'][:,i])
    vectUInvProjAtIntegPoints[i] = phiAtIntegPoint.dot(unstructuredMeshClipped.nodeFields['UInvProj'][:,i])

normDelta = np.einsum('ij,ij,j', vectDeltaAtIntegPoints, vectDeltaAtIntegPoints, integrationWeights, optimize = True)
normUexact = np.einsum('ij,ij,j', vectUexactAtIntegPoints, vectUexactAtIntegPoints, integrationWeights, optimize = True)
normUInvProj = np.einsum('ij,ij,j', vectUInvProjAtIntegPoints, vectUInvProjAtIntegPoints, integrationWeights, optimize = True)

print("normDelta =", normDelta)
print("normUexact =", normUexact)
print("normUInvProj =", normUInvProj)

print("rel normDelta =", normDelta/normUexact)


unstructuredMeshClipped.nodeFields['deltaRel'] = (unstructuredMeshClipped.nodeFields['Uexact']-unstructuredMeshClipped.nodeFields['UInvProj'])/normUexact


from BasicTools.IO import XdmfWriter as XW
#write solution in xdmf
writer = XW.XdmfWriter('ClippedMesh.xmf')
writer.SetTemporal(True)
writer.SetHdf5(False)
writer.Open()
writer.Write(unstructuredMeshClipped,PointFields=[unstructuredMeshClipped.nodeFields['Uexact'][:,0],unstructuredMeshClipped.nodeFields['Uexact'][:,1],unstructuredMeshClipped.nodeFields['UInvProj'][:,0],unstructuredMeshClipped.nodeFields['UInvProj'][:,1],unstructuredMeshClipped.nodeFields['delta'][:,0],unstructuredMeshClipped.nodeFields['delta'][:,1],unstructuredMeshClipped.nodeFields['deltaRel'][:,0],unstructuredMeshClipped.nodeFields['deltaRel'][:,1]], PointFieldsNames=['UexactX','UexactY','UInvProjX','UInvProjY','deltaX','deltaY','deltaRelX','deltaRelY'])
writer.Close()























