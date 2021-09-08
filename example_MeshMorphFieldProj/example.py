#from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh
from BasicTools.Containers.UnstructuredMeshModificationTools import Morphing, ComputeSkin

from BasicTools.IO import XdmfWriter as XW
import numpy as np

structuredMesh = ConstantRectilinearMesh(dim=3)
structuredMesh.SetDimensions([8,12,16])
structuredMesh.SetOrigin([0.,0.,0.])
structuredMesh.SetSpacing([0.1,0.1,0.1])

unstructuredMesh = CreateMeshFromConstantRectilinearMesh(structuredMesh,ofTetras=True)

print(unstructuredMesh)

ComputeSkin(unstructuredMesh, inplace=True)

print(unstructuredMesh)

ids = unstructuredMesh.GetElementsInTag("Skin")

print(ids)

from BasicTools.Containers.Filters import NodeFilter, ElementFilter 

nodeFilter = NodeFilter(unstructuredMesh)
nodeFilter.AddETag("Skin") #696 alone
nodeFilter.AddZone(lambda p: 1-p[:,2]) #864 alone   #both: 384

idsToTreat = nodeFilter.GetIdsToTreat()
print("len idsToTreat =", len(idsToTreat))

nodalTag = unstructuredMesh.GetNodalTag("ImposedDisp")
nodalTag.AddToTag(idsToTreat)


posIdsToTreat = unstructuredMesh.GetPosOfNodes()[idsToTreat]

dispIdsToTreat = np.zeros(posIdsToTreat.shape)
dispIdsToTreat[:,0] = (posIdsToTreat[:,2] - 1.) * (posIdsToTreat[:,0] - 1.)
dispIdsToTreat[:,1] = (posIdsToTreat[:,2] - 1.) * (posIdsToTreat[:,1] - 1.)

print("pos =", posIdsToTreat)
print("pos =", dispIdsToTreat)

writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='init.xdmf')
writer.Write(unstructuredMesh, PointFields = [], PointFieldsNames = [])
writer.Close()

unstructuredMesh.nodes[idsToTreat] += dispIdsToTreat

writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='onlySkinMorphed.xdmf')
writer.Write(unstructuredMesh, PointFields = [], PointFieldsNames = [])
writer.Close()

unstructuredMesh.nodes[idsToTreat] -= dispIdsToTreat

unstructuredMesh.SetNodes(Morphing(unstructuredMesh, dispIdsToTreat, idsToTreat))
print(unstructuredMesh)





unstructuredMesh.ComputeBoundingBox()

structuredMeshBB = ConstantRectilinearMesh(dim=3)
n = 10
spacing = [(unstructuredMesh.boundingMax[i]-unstructuredMesh.boundingMin[i])/(n-1) for i in range(3)]
structuredMeshBB.SetDimensions([n,n,n])
structuredMeshBB.SetOrigin(unstructuredMesh.boundingMin)
structuredMeshBB.SetSpacing(spacing)


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


writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='boundingBox.xdmf')
writer.Write(structuredMeshBB, PointFields = [testField], PointFieldsNames = ["testField"])
writer.Close()





from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh
from BasicTools.FE.FETools import PrepareFEComputation
from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
from BasicTools.FE.Fields.FEField import FEField


unstructuredMeshBB = CreateMeshFromConstantRectilinearMesh(structuredMeshBB)


print("unstructuredMeshBB =", unstructuredMeshBB)

space, numberings, _offset, _NGauss = PrepareFEComputation(unstructuredMeshBB, numberOfComponents = 1)
inputFEField = FEField(name="testField",mesh=unstructuredMeshBB,space=space,numbering=numberings[0])
operator, status = GetFieldTransferOp(inputFEField, unstructuredMesh.nodes, method = "Interp/Clamp", verbose=True)


print(testField.shape)
print(operator.shape)
projectedTestField = operator.dot(testField)


writer = XW.XdmfWriter()
writer.SetXmlSizeLimit(0)
writer.SetBinary(True)
writer.Open(filename='morphed.xdmf')
writer.Write(unstructuredMesh, PointFields = [projectedTestField], PointFieldsNames = ["projectedTestField"])
writer.Close()






from BasicTools.FE.Integration import IntegrateGeneral
from BasicTools.FE.SymWeakForm import GetField
from BasicTools.FE.SymWeakForm import GetTestField

F = GetField("F",1)
Tt = GetTestField("T",1)

wf = F.T*Tt


space, numberings, _offset, _NGauss = PrepareFEComputation(unstructuredMesh, numberOfComponents = 1)


field = FEField("F",mesh=unstructuredMesh,space=space,numbering=numberings[0], data = projectedTestField)
unkownField = FEField("T",mesh=unstructuredMesh,space=space,numbering=numberings[0])

elFilter = ElementFilter(unstructuredMesh)
elFilter.AddTag("Skin")

K, F = IntegrateGeneral(mesh=unstructuredMesh,
                    wform=wf,
                    constants={},
                    fields=[field],
                    unkownFields=[unkownField],
                    integrationRuleName="LagrangeP1",
                    elementFilter=elFilter)


#print(K)
print(K.shape)
print(K.dot(projectedTestField))
print(np.linalg.norm(F))
print(np.sum(F))





#print("pos =", unstructuredMesh.GetPosOfNodes()[idsToTreat])

""" ##### method for computing the deform mesh knowing displacement of some nodes #######################################
 ##### BC is the known displacement in a numpy array (shape [number_of_of_known_nodes,3]) ############################
 ##### bord_tot contains the ids of known nodes (list of ids or boolean array) 


    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)
    BC = np.empty((mesh.GetNumberOfNodes(),3),dtype=float)
    bord_tot = np.empty(mesh.GetNumberOfNodes(),dtype=int)
    cpt = 0
    print(mesh)
    for name,data in mesh.elements.items():

        if ElementNames.dimension[name] != 2:
            continue

        ids = data.GetNodesIdFor(data.GetTag("X0").GetIds())
        print(ids)
        bord_tot[cpt:cpt+len(ids)] = ids
        BC[cpt:cpt+len(ids),:] = 0
        cpt += len(ids)

        ids = data.GetNodesIdFor(data.GetTag("X1").GetIds())
        print(ids)
        bord_tot[cpt:cpt+len(ids)] = ids
        BC[cpt:cpt+len(ids),:] = [[0,0,10]]
        cpt += len(ids)

    BC = BC[0:cpt,:]
    bord_tot = bord_tot[0:cpt]

    new_p1 = Morphing(mesh, BC,bord_tot)
    new_p2 = Morphing(mesh, BC,bord_tot,rayon= 20. )

    new_p0 = np.copy(mesh.nodes)
    new_p0[bord_tot,:] += BC
    mesh.nodeFields["morph0"] = new_p0
    mesh.nodeFields["morph1"] = new_p1
    mesh.nodeFields["morph2"] = new_p2










elFilter = ElementFilter(unstructuredMesh)
elFilter.AddTag("Skin")

nodeIDs = set()
for name,data,ids in elFilter:
    print("name,data,ids =", name,data.connectivity,ids)
    for nodeIDsEl in data.connectivity:
        for nodeID in nodeIDsEl:
            nodeIDs.add(nodeID)

unstructuredMesh.nodesTags.CreateTag('Skin')
print("unstructuredMesh.nodesTags =", unstructuredMesh.nodesTags)
print("unstructuredMesh.nodesTags['Skin'] =", unstructuredMesh.nodesTags['Skin'])

unstructuredMesh.nodesTags['Skin'].AddToTag(list(nodeIDs))

print("unstructuredMesh.nodesTags['Skin'].GetIds() =", unstructuredMesh.nodesTags['Skin'].GetIds())

print("unstructuredMesh =", unstructuredMesh)


print("noFilter =", noFilter)


print("noFilter.GetIdsToTreat() =", noFilter.GetIdsToTreat())

#elFilter.AddETag(Skin)



"""
