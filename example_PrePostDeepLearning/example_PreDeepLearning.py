from BasicTools.IO import XdmfReader as XR
import numpy as np

reader = XR.XdmfReader(filename = "jet_flame.sol000829.xmf")
reader.Read();
dom = reader.xdmf.GetDomain(0);

grid = dom.GetGrid(0)
mesh = grid.GetSupport()

PointFieldsNames = grid.GetPointFieldsNames()
#allFields = np.array(grid.GetPointFields())
U = grid.GetPointFields()[0][:,0:2]


print("mesh =", mesh)
print("PointFieldsNames =", PointFieldsNames)
print("U =", U)
print("U.shape =", U.shape)

Nx = 48
Ny = 48

originMesh = [0., np.min(mesh.nodes[:,1])]
coordMax = [np.max(mesh.nodes[:,0]), np.max(mesh.nodes[:,1])]

Lx = coordMax[0] - originMesh[0]  
Ly = coordMax[1] - originMesh[1]  

print(Lx,Ly)
    
#Create Projected mesh  
from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh    
rectMesh = ConstantRectilinearMesh(dim=2)
rectMesh.SetDimensions([Nx,Ny])
rectMesh.SetSpacing([Lx/(Nx-1), Ly/(Ny-1)])
rectMesh.SetOrigin([originMesh[0],originMesh[1]])

print("rectMesh =", rectMesh)



import BasicTools.Containers.UnstructuredMeshCreationTools as UMCT
from BasicTools.FE.FETools import PrepareFEComputation
from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
from BasicTools.FE.Fields.FEField import FEField

rectMesh2 = UMCT.CreateMeshFromConstantRectilinearMesh(rectMesh)
space, numberings, _offset, _NGauss = PrepareFEComputation(mesh, numberOfComponents = 1)
inputFEField = FEField(name="toto",mesh=mesh,space=space,numbering=numberings[0])
methods = ["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"]
method = methods[2]
operator, status = GetFieldTransferOp(inputFEField, rectMesh2.nodes, method = method, verbose=True)

nbeGrids = 1
res = np.empty((nbeGrids,2,Nx,Ny))
for i in range(nbeGrids):
        projectedField = operator.dot(U[:,0])
        res[i,0,:,:] = projectedField.reshape((Nx,Ny))
        projectedField = operator.dot(U[:,1])
        res[i,1,:,:] = projectedField.reshape((Nx,Ny))

print("res =", res)


from BasicTools.IO import XdmfWriter as XW
writer = XW.XdmfWriter('TestProj.xmf')
writer.SetTemporal(True)
writer.SetHdf5(True)
writer.Open()
for i in range(nbeGrids):
    writer.Write(rectMesh,PointFields=[res[i,0], res[i,1]], PointFieldsNames=["UX","UY"]);
writer.Close()




