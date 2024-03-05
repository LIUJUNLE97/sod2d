#!/usr/bin/env python
#
# Ansa2sod
#
# Export a mesh from Ansa to Sod2D format.
#
# Last rev: 27/10/2023

import h5py, numpy as np

import ansa
from ansa import base, constants, mesh, utils

deck = constants.CGNS
###############################
#####   Case Definition   #####
###############################
path = "/home/sgome1/HAMR/TEST/"
basename = 'tgv4'
scaleFactor = 1.0
perDist = np.pi*2
###############################
##### End Case Definition #####
###############################

class ordering(object):

    def __init__(self,pOrder,fmt='GMSH'):
        '''
        Generates the ordering of the element and its faces given a format
        '''
        self._pOrder = pOrder
        self._nedge  = self._pOrder+1
        self._nnode  = self._nedge**3

        if fmt =='GMSH':
            self._quad_order_edges = np.array([0,1,1,2,2,3,3,0],dtype='int32').reshape((4,2))
            self._hex_order_edges = np.array([0, 1, 0, 3, 0, 4, 1, 2, 1, 5, 2, 3,
            2, 6, 3, 7, 4, 5, 4, 7, 5, 6, 6, 7],dtype='int32').reshape((12,2))
            self._hex_order_faces = np.array([0, 3, 2, 1, 0, 1, 5, 4, 0, 4, 7, 3,
             1, 2, 6, 5, 2, 3, 7, 6, 4, 5, 6, 7],dtype='int32').reshape((6,4))
            self._indexTable = self.hexaHOtable(self._pOrder)
            self.ijkTable()
        elif fmt == 'CGNS':
            self._quad_order_edges = np.array([0,1,1,2,2,3,3,0],dtype='int32').reshape((4,2))
            self._hex_order_edges = np.array([0,1,1,2,2,3,3,0,0,4,1,5,2,6,3,7,4,5,5,6,6,7,7,4],dtype='int32').reshape((12,2))
            # self._hex_order_faces = np.array([0,1,2,3,0,1,5,4,1,2,6,5,2,3,7,6,0,4,7,3,4,5,6,7],dtype='int32').reshape((6,4))
            self._hex_order_faces = np.array([0,1,2,3,0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7,4,5,6,7],dtype='int32').reshape((6,4))
            self._indexTable = self.cgnsHexaTable(self._pOrder)
            self.ijkTable()
            self.nodesAnsaIJK()
        else:
            print("This format is not available")
            return

    def hexaHOtable(self,p):
        indexTable = np.zeros(((p+1)**3,3),dtype='int32')
        if p > 0:
            #                  i, j, k
            indexTable[1,:] = [p, 0, 0]
            indexTable[2,:] = [p, p, 0]
            indexTable[3,:] = [0, p, 0]
            indexTable[4,:] = [0, 0, p]
            indexTable[5,:] = [p, 0, p]
            indexTable[6,:] = [p, p, p]
            indexTable[7,:] = [0, p, p]
        if p > 1:
            # Generate high-order edges
            inode = 8
            for iedge in range(len(self._hex_order_edges)):
                i0 = self._hex_order_edges[iedge,0]
                i1 = self._hex_order_edges[iedge,1]
                u = (indexTable[i1,:] - indexTable[i0,:])/p
                for i in range(1,p):
                    indexTable[inode,:] = indexTable[i0,:] + i*u
                    inode += 1
            # Generate a generic high-order face with p = p-2
            tableFace = self.quadHOtable(p-2)
            tableFace+=1
            # Generate faces interior nodes
            for iface in range(len(self._hex_order_faces)):
                i0 = self._hex_order_faces[iface,0]
                i1 = self._hex_order_faces[iface,1]
                i3 = self._hex_order_faces[iface,3]
                u = (indexTable[i1,:] - indexTable[i0,:])/p
                v = (indexTable[i3,:] - indexTable[i0,:])/p
                for i in range((p-1)**2):
                    indexTable[inode,:] = indexTable[i0,:] + u*tableFace[i,0] + v*tableFace[i,1]
                    inode += 1
            # Generate volume nodes
            tableVolume = self.hexaHOtable(p-2)
            tableVolume+=1
            indexTable = self.joinTables(inode,tableVolume,indexTable)
        return indexTable

    def quadHOtable(self,p):
        indexTable = np.zeros(((p+1)**2,2),dtype='int32')
        tableFace  = np.zeros(((p-1)**2,2),dtype='int32')
        indexTable[0,:] = [0,0]
        if p > 0:
            indexTable[1,:] = [p,0]
            indexTable[2,:] = [p,p]
            indexTable[3,:] = [0,p]
        if p > 1:
            inode = 4
            for iedge in range(len(self._quad_order_edges)):
                i0 = self._quad_order_edges[iedge,0]
                i1 = self._quad_order_edges[iedge,1]
                u = (indexTable[i1,:] - indexTable[i0,:])/p
                for i in range(1,p):
                    indexTable[inode,:] = indexTable[i0,:] + i*u
                    inode = inode + 1
            # breakpoint()
            tableFace = self.quadHOtable(p-2)
            tableFace+=1
            indexTable = self.joinTables(inode,tableFace,indexTable)
        return indexTable

    def joinTables(self,indexDesti,table1,table2):
        j = indexDesti
        for i in range(len(table1)):
            table2[j,:] = table1[i,:]
            j+=1
        return table2

    def cgnsHexaTable(self,p):
        indexTable = np.zeros(((p+1)**3,3),dtype='int32')
        if p > 0:
            #                  i, j, k
            indexTable[1,:] = [p, 0, 0]
            indexTable[2,:] = [p, p, 0]
            indexTable[3,:] = [0, p, 0]
            indexTable[4,:] = [0, 0, p]
            indexTable[5,:] = [p, 0, p]
            indexTable[6,:] = [p, p, p]
            indexTable[7,:] = [0, p, p]
        if p > 1:
            # Generate high-order edges
            inode = 8
            for iedge in range(len(self._hex_order_edges)):
                i0 = self._hex_order_edges[iedge,0]
                i1 = self._hex_order_edges[iedge,1]
                u = (indexTable[i1,:] - indexTable[i0,:])/p
                for i in range(1,p):
                    indexTable[inode,:] = indexTable[i0,:] + i*u
                    inode += 1
            # Generate a generic high-order face with p = p-2
            tableFace = self.cgnsQuadTable(p-2)
            tableFace+=1
            # Generate faces interior nodes
            for iface in range(len(self._hex_order_faces)):
                i0 = self._hex_order_faces[iface,0]
                i1 = self._hex_order_faces[iface,1]
                i3 = self._hex_order_faces[iface,3]
                u = (indexTable[i1,:] - indexTable[i0,:])/p
                v = (indexTable[i3,:] - indexTable[i0,:])/p
                for i in range((p-1)**2):
                    indexTable[inode,:] = indexTable[i0,:] + u*tableFace[i,0] + v*tableFace[i,1]
                    inode += 1
            # Generate volume interior nodes
            i0 = self._hex_order_faces[0,0]
            i1 = self._hex_order_faces[0,1]
            i3 = self._hex_order_faces[0,3]
            u = (indexTable[i1,:] - indexTable[i0,:])/p
            v = (indexTable[i3,:] - indexTable[i0,:])/p
            for iface in range(1,p):
                for i in range((p-1)**2):
                    indexTable[inode,:] = indexTable[i0,:] + u*tableFace[i,0] + v*tableFace[i,1]
                    indexTable[inode,2] = iface
                    inode += 1
        return indexTable

    def cgnsQuadTable(self,p):
        quad_order_edges = np.array([0,1,1,2,2,3,3,0],dtype='int32').reshape((4,2))
        indexCorners = np.array([[0,0],[p,0],[p,p],[0,p]])
        indexTable = np.zeros(((p+1)**2,2),dtype='int32')
        indexTable[0,:] = [0,0]
        if p > 0:
            inode = 1
            for iedge in range(3):
                i0 = self._quad_order_edges[iedge,0]
                i1 = self._quad_order_edges[iedge,1]
                u = (indexCorners[i1,:] - indexCorners[i0,:])/p
                for i in range(1,p+1):
                    indexTable[inode,:] = indexCorners[i0,:] + i*u
                    inode = inode + 1
            p_ = p-1
            while p_ > 0:
                for j in range(2):
                    iedge+=1
                    edge = iedge%4
                    i0 = self._quad_order_edges[edge,0]
                    i1 = self._quad_order_edges[edge,1]
                    u = (indexCorners[i1,:] - indexCorners[i0,:])/p
                    for i in range(1,p_+1):
                        indexTable[inode,:] = indexTable[inode-1,:] + u
                        inode = inode + 1
                p_-=1
        return indexTable

    def ijkTable(self):
        self._ijkTable = np.zeros((self._nedge,self._nedge,self._nedge),dtype='int32')
        for pIndex in range(self._nnode):
            i = self._indexTable[pIndex,0]
            j = self._indexTable[pIndex,1]
            k = self._indexTable[pIndex,2]
            self._ijkTable[i,j,k] = pIndex

    def get_2ijk(self):
        order = np.zeros((self._nnode),dtype='int32')
        for k in range(self._nedge):
            for j in range(self._nedge):
                for i in range(self._nedge):
                    inode = k*self._nedge**2 + j*self._nedge + i
                    order[inode] = self._ijkTable[i,j,k]
        return order

    def getIJ2_(self):
        ijOrder = self.quadHOtable(self._pOrder)
        order = ijOrder[:,0]+ijOrder[:,1]*self._nedge
        return order

    def getIJK2_(self):
        order = np.zeros((self._nnode),dtype='int32')
        for inode in range(self._nnode):
            i = self._indexTable[inode,0]
            j = self._indexTable[inode,1]
            k = self._indexTable[inode,2]
            order[inode] = k*self._nedge**2 + j*self._nedge + i
        return order

    def nodesAnsa(self, nodes):
        nodesTags = []
        for node in nodes:
            nodesTags.append('N{}'.format(node))
        return nodesTags
    
    def nodesAnsaIJK(self):
        nodes = np.arange(1,self._nnode+1,dtype='int32')
        order = self.get_2ijk()
        self._nodesIJK = nodes[order]

    def facesIndecesREF(self):
        self._facesIndecesREF = np.zeros((6,self._nedge**2),dtype='int32') # Left,Front,Bottom,Rigth,Back,Top
        ijk = np.zeros((3,),dtype='int32')
        axises = np.arange(3)
        values = np.array([0,self._pOrder])
        face = 0
        for value1 in values:
            for axis1 in axises:
                for axis2 in axises:
                    if axis2 != axis1:
                        for axis3 in axises:
                            if axis3 != axis2 and axis3 != axis1 and axis3 > axis2:
                                ID = 0
                                for value2 in range(self._nedge):
                                    for value3 in range(self._nedge):
                                        ijk[axis1]=value1
                                        ijk[axis2]=value2
                                        ijk[axis3]=value3
                                        self._facesIndecesREF[face,ID] = ijk[0]+ijk[1]*self._nedge+ijk[2]*self._nedge**2
                                        ID+=1
                face += 1

    def facesIndeces(self):
        self.facesIndecesREF()
        self._facesIndeces = np.zeros((6,self._nedge**2),dtype='int32') # Left,Front,Bottom,Rigth,Back,Top
        for iface,face in enumerate(self._facesIndecesREF):
            self._facesIndeces[iface,:] = self._nodesIJK[face]
        return self._facesIndeces

def getPolynomialOrder(info):
    hexaOrders = dict(HEXA = 1,HEXA_27 = 2,HEXA_64 = 3,HEXA_125 = 4,)
    hexaType = list(info['ELEMENT'].children['SOLID'].children.keys())[0]
    pOrder = hexaOrders[hexaType]
    return pOrder

# Collect only shell PIDs in which "USE_IN_MODEL" option is YES, i.e. for avoid collecting top_cap PID.
def checkShells(base, deck):
    shells = []
    pershells = []
    usedPshells = []
    periodicPshells = []
    pshells = base.CollectEntities(deck, None, 'SHELL_PROPERTY')
    for pshell in pshells:
        vals = base.GetEntityCardValues(deck, pshell, ['USE_IN_MODEL', 'TYPE'])
        if vals['USE_IN_MODEL'] == 'YES':
            usedPshells.append(pshell)
            pidShells = base.CollectEntities(deck, pshell, 'SHELL', recursive = True)
            shells = shells + pidShells 
        elif vals['TYPE'] != 'Interior':
            periodicPshells.append(pshell)
            pidShells = base.CollectEntities(deck, pshell, 'SHELL', recursive = True)
            pershells = pershells + pidShells
    return(len(shells), usedPshells, shells, len(pershells), periodicPshells, pershells)

# Renumber nodes
def renumberNodes(base, deck, nnodes):
    i=1
    j=1
    while j <= nnodes:
        try:
            node = base.GetEntity(deck, 'NODE', i)
            base.SetEntityId(node, j, True, False)
            j+=1
        except:
            pass
        i+=1

# Write nodes
def writeNodes(base, deck, h5file, scaleFactor, nnodes):
    xyz = np.zeros((nnodes,3))
    for i in range(1, nnodes+1):
        node = base.GetEntity(deck, 'NODE', i)
        vals = base.GetEntityCardValues(deck, node, ['X', 'Y', 'Z'])
        xyz[i-1,:] = np.array([vals['X'],vals['Y'],vals['Z']])*scaleFactor
    nodes_dset = h5file.create_dataset('coords',(nnodes,3),dtype='f8',data=xyz,chunks=True,maxshape=(nnodes,3))
    del xyz

# Write boundary elements
def writeBounds(base, deck, h5file, nshells, shells, pOrder, nameDataSet):
    NOD_QUAD = ['N4','N3','N2','N1']
    # NOD_QUAD = ['N4','N3','N2','N1','N5','N6','N7','N8','N9']
    nnodeB = len(NOD_QUAD)
    bounds = np.zeros((nshells,nnodeB))
    for i,shell in enumerate(shells):
        vals = base.GetEntityCardValues(deck, shell, ['type'])
        type = vals['type']
        if type == 'QUAD':
            qvals = base.GetEntityCardValues(deck, shell, NOD_QUAD)
            for j in range(len(NOD_QUAD)):
                bounds[i,j] = qvals[NOD_QUAD[j]]
    bounds_dset = h5file.create_dataset(nameDataSet,(nshells,nnodeB),dtype='i8',data=bounds,
        chunks=True,maxshape=(nshells,nnodeB))
    del bounds

def getElemFaceNodes(deck, base, elemID, shell, pOrder):
    NOD_QUAD = ['N4','N3','N2','N1']
    faceElemCorners = [
                    ['N4','N1','N5','N8'], # Left
                    ['N1','N5','N2','N6'], # Front
                    ['N1','N2','N3','N4'], # Bottom
                    ['N2','N3','N7','N6'], # Right
                    ['N3','N4','N8','N7'], # Back
                    ['N5','N6','N7','N8']] # Top
    faceElemNodes = ordering(pOrder,'CGNS').facesIndeces()
    bound = np.zeros(((pOrder+1)**2,),dtype='int32')
    elem = base.GetEntity(deck,'SOLID', elemID)
    cornerEntities = base.CollectEntities(deck, shell, 'NODE', recursive = True)
    cornerNodes = np.zeros((len(NOD_QUAD)), dtype='int')
    for i,node in enumerate(cornerEntities):
        cornerNodes[i] = node._id
    cornerNodes = np.sort(cornerNodes)
    for i, face in enumerate(faceElemCorners):
        faceCorners = base.GetEntityCardValues(deck,elem,faceElemCorners[i])
        elemFaceCornerNodes = np.zeros((len(NOD_QUAD)), dtype='int')
        for j,node in enumerate(faceCorners):
            elemFaceCornerNodes[j] = faceCorners[face[j]]
        elemFaceCornerNodes = np.sort(elemFaceCornerNodes)
        if np.all(elemFaceCornerNodes == cornerNodes):
            faceID = i
            faceNodes = ordering(pOrder,'CGNS').nodesAnsa(faceElemNodes[faceID])
            vals = base.GetEntityCardValues(deck,elem,faceNodes)
            for j in range(len(faceNodes)):
                bound[j] = vals[faceNodes[j]]
            break
    return bound

def writeBoundsAlt(base, deck, h5file, nshells, shells, pOrder, nameDataSet):
    nnodeB = (pOrder+1)**2
    bounds = np.zeros((nshells,nnodeB))
    if nshells > 0:
        elems = base.CollectEntities(deck, None, 'SOLID', recursive = True)
        pairs = mesh.MatchShellsAndSolids(shells, elems)
        del elems
        d={}
        for i in range(int(len(pairs)/2)):
            d[pairs[i*2-2]._id] = pairs[i*2-1]._id
        for iShell,shell in enumerate(shells):
            bound = getElemFaceNodes(deck, base, d[shell._id], shell, pOrder)
            bounds[iShell,:] = bound
        order = ordering(pOrder,'GMSH').getIJ2_()
        bounds = bounds[:,order]
    bounds_dset = h5file.create_dataset(nameDataSet,(nshells,nnodeB),dtype='i8',data=bounds,
    chunks=True,maxshape=(nshells,nnodeB))
    

    if nameDataSet == 'periodicFaces':
        perNodes = bounds.flatten().astype('int32')
        del bounds
        return perNodes
    else:
        del bounds

# Write boundary codes
def writeBC(base, deck, h5file, usedPshells, nbouns):
    bound_codes = np.zeros((nbouns,))
    iShell = 0
    for pshell in usedPshells:
        code = base.GetEntityCardValues(deck, pshell, ['ZONE_ID'])['ZONE_ID']
        shells = base.CollectEntities(deck, pshell, 'SHELL', recursive = True)
        for shell in shells:
            bound_codes[iShell] = code
            iShell+=1
    boundId_dset = h5file.create_dataset('boundFacesId',(nbouns,),dtype='i8',data=bound_codes,
        chunks=True,maxshape=(nbouns))
    del bound_codes

# Write pair of periodic nodes
def periodicPair(base, deck, h5file, periodicPshells, nbouns, dims_group):
    parents = []
    childs = []
    npairs = 0
    for pshell in periodicPshells:
        faces = base.CollectEntities(deck, pshell,'FACE', recursive = True)
        for face in faces:
            vals = base.GetEntityCardValues(deck, face, ['ID','Child Link Faces'])
            parentID = vals['ID']
            childID = vals['Child Link Faces']
            if len(childID) > 0:
                parents.append(parentID)
                childs.append(int(childID))
                npairs += len(base.CollectEntities(deck, face,'NODE', recursive = True))
    
    dset = dims_group.create_dataset('numPeriodicLinks',(1,),dtype='i8',data=npairs)
    pairs = np.zeros((npairs,2))
    ipair = 0
    for parentID, childID in zip(parents,childs):
        parentFace = base.GetEntity(deck, 'FACE', parentID)
        childFace = base.GetEntity(deck, 'FACE', childID)
        parentNodes = base.CollectEntities(deck, parentFace,'NODE', recursive = True)
        childNodes = base.CollectEntities(deck, childFace,'NODE', recursive = True)[::-1]
        for parent, child in zip(parentNodes,childNodes):
            pairs[ipair,:] = [parent._id,child._id]
            ipair+=1
    pairs_dset = h5file.create_dataset('periodicLinks',(npairs,2),dtype='i8',data=pairs,
        chunks=True,maxshape=(npairs,2))
    del pairs

def periodicPairAlt(base, deck, h5file, periodicPshells, dims_group, perDist,pOrder):
    parentPIDS = periodicPshells[:len(periodicPshells)//2]
    childPIDS = periodicPshells[len(periodicPshells)//2:]
    for parentPID, childPID in zip(parentPIDS,childPIDS):
        parentShells = base.CollectEntities(deck, parentPID, 'SHELL', recursive = True)
        childShells = base.CollectEntities(deck, childPID, 'SHELL', recursive = True)
        elems = base.CollectEntities(deck, None, 'SOLID', recursive = True)
        parentPairs = mesh.MatchShellsAndSolids(parentShells, elems)
        childPairs = mesh.MatchShellsAndSolids(childShells, elems)
        del elems
        parentDic={}
        childDic={}
        for i in range(int(len(parentPairs)/2)):
            parentDic[parentPairs[i*2-2]._id] = parentPairs[i*2-1]._id
            childDic[childPairs[i*2-2]._id] = childPairs[i*2-1]._id
        parentNodes = np.zeros((len(parentShells),(pOrder+1)**2),dtype='int32')
        childNodes = np.zeros((len(childShells),(pOrder+1)**2),dtype='int32')
        i = 0
        for parentShell,childShell in zip(parentShells,childShells):
            parentShellNodes = getElemFaceNodes(deck, base, parentDic[parentShell._id], parentShell, pOrder)
            childShellNodes = getElemFaceNodes(deck, base, childDic[childShell._id], childShell, pOrder)
            parentNodes[i,:] = parentShellNodes
            childNodes[i,:] = childShellNodes
            i+=1
        uniqueParentNodes = np.unique(parentNodes.flatten())
        uniqueChildNodes = np.unique(childNodes.flatten())

        xyzParent = np.zeros((len(uniqueParentNodes),3))
        xyzChild = np.zeros((len(uniqueChildNodes),3))
        for i,idNode in enumerate(uniqueParentNodes):
            node = base.GetEntity(deck, 'NODE', idNode)
            vals = base.GetEntityCardValues(deck, node, ['X', 'Y', 'Z'])
            xyzParent[i,:] = [vals['X'],vals['Y'],vals['Z']]
        for i,idNode in enumerate(uniqueChildNodes):
            node = base.GetEntity(deck, 'NODE', idNode)
            vals = base.GetEntityCardValues(deck, node, ['X', 'Y', 'Z'])
            xyzChild[i,:] = [vals['X'],vals['Y'],vals['Z']]

        localPairs = np.zeros((len(uniqueParentNodes),2),dtype = 'int32')
        for i,xyz in enumerate(xyzParent):
            dist = np.abs(np.sqrt(np.sum((xyzChild-xyz)**2,axis=1))-perDist)
            iChild = uniqueChildNodes[np.argmin(dist)]
            iParent = uniqueParentNodes[i]
            localPairs[i,:] = np.array([iParent,iChild])
            xyzChild = np.delete(xyzChild,np.argmin(dist),0)
            uniqueChildNodes = np.delete(uniqueChildNodes,np.argmin(dist))
        if 'pairs' in locals():
            pairs = np.vstack((pairs,localPairs))
        else:
            pairs = localPairs

    dset = dims_group.create_dataset('numPeriodicLinks',(1,),dtype='i8',data=len(pairs))
    pairs_dset = h5file.create_dataset('periodicLinks',(len(pairs),2),dtype='i8',data=pairs,
        chunks=True,maxshape=(len(pairs),2))
    del pairs

# Write solids
def writeElems(base, deck, h5file, nelems,pOrder):
    nodes = ordering(pOrder,'CGNS')._nodesIJK
    NOD_HEXA = ordering(pOrder,'CGNS').nodesAnsa(nodes)
    nnodeE = len(NOD_HEXA)
    connec = np.zeros((nelems,nnodeE))
    elems = base.CollectEntities(deck, None, 'SOLID', recursive = True)
    for i, elem in enumerate(elems):
        hvals = base.GetEntityCardValues(deck, elem, NOD_HEXA)
        for j in range(len(NOD_HEXA)):
            connec[i,j] = hvals[NOD_HEXA[j]]
    order = ordering(pOrder,'GMSH').getIJK2_()
    connec = connec[:,order]
    connec_dset = h5file.create_dataset('connec',(nelems,nnodeE),dtype='i8',data=connec,
        chunks=True,maxshape=(nelems,nnodeE))
    del elems

# Write important information about the domain, geometry and boundary codes files
def writeInfo(base, deck, scaleFactor, usedPshells, path):
    o = open(path+'info.dat', 'w')
    o.write('Ansa2sod2d' + '\n\n')
    o.write('This output was generated from the following file ' + '\n' + base.DataBaseName() + ':\n\n')
    o.write('SCALE FACTOR = ' + str(scaleFactor) + '\n\n')
    if len(usedPshells)>0:
        o.write('BOUNDARY_CODES'+'\n')
        writeBoundaryName(base, deck, o, usedPshells)
        o.write('END_BOUNDARY_CODES'+'\n')
        o.close()

# Write the name of the differents boundaries and their codes
def writeBoundaryName(base, deck, out, usedPshells):
    i=0
    for pshell in usedPshells:
        vals = base.GetEntityCardValues(deck, pshell, ['Name','__id__'])
        pshellname = vals['Name']
        pshellID = vals['__id__']
        shellLine = '{:4d} {}\n'.format(pshellID,pshellname)
        out.write(shellLine)

def main():
    # Open HDF5 file
    h5filename = path+basename+'.h5'
    print(h5filename)
    h5file = h5py.File(h5filename,'w')
    dims_group = h5file.create_group('dims')
    # Compresses unused F.E. entities and deleted geometrical entities.
    base.Compress('')
    # Get general mesh info
    info = utils.DatabaseBrowserInfo(deck)
    nelems = info['ELEMENT'].children['SOLID'].total
    nnodes = info['NODE'].total
    pOrder = getPolynomialOrder(info)
    nbouns, usedPshells, shells, nper, periodicPshells, pershells = checkShells(base, deck)
    # Write dims
    dset = dims_group.create_dataset('order',(1,),dtype='i8',data=pOrder)
    dset = dims_group.create_dataset('numNodes',(1,),dtype='i8',data=nnodes)
    elems_dset = dims_group.create_dataset('numElements',(1,),dtype='i8',data=nelems)
    bound_dset = dims_group.create_dataset('numBoundaryFaces',(1,),dtype='i8',data=nbouns)
    per_dset   = dims_group.create_dataset('numPeriodicFaces',(1,),dtype='i8',data=nper)
    # Renumbering nodes
    renumberNodes(base, deck, nnodes)
    # Writing nodes coordinates
    writeNodes(base, deck, h5file, scaleFactor, nnodes)
    # Writing boundaries
    # writeBounds(base, deck, h5file, nbouns, shells, pOrder, 'boundFaces')
    writeBoundsAlt(base, deck, h5file, nbouns, shells, pOrder, 'boundFaces')
    # Writing boundary codes
    writeBC(base, deck, h5file, usedPshells, nbouns)
    # Writing periodic boundaries
    if len(periodicPshells) > 0:
        # writeBounds(base, deck, h5file, nper, pershells, pOrder, 'periodicFaces')
        perNodes = writeBoundsAlt(base, deck, h5file, nper, pershells, pOrder, 'periodicFaces')
        # Writing periodic pairs
        # periodicPair(base, deck, h5file, periodicPshells, nper, dims_group)
        periodicPairAlt(base, deck, h5file, periodicPshells, dims_group, perDist, pOrder)

    # Writing elements
    writeElems(base, deck, h5file, nelems, pOrder)

    
    print('Writing .info.dat...' + '\n')
    writeInfo(base, deck, scaleFactor, usedPshells, path)
    print('.info.dat written...' + '\n')

    print("Done!!!")


if __name__ == '__main__':
    main()