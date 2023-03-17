#!/usr/bin/env python
#
# gmsh2sod2d
#
# Export a mesh from GMSH to Sod2D format.
#
# Last rev: 14/03/2023
from __future__ import print_function, division

# Please do not delete this part otherwise it will not work
# you have been warned after a long weekend of debugging
import mpi4py
mpi4py.rc.recv_mprobe = False
from mpi4py import MPI

#import os, sys, itertools, argparse, numpy as np
import h5py, argparse, numpy as np
#import #pyAlya

MPI_RANK = MPI.COMM_WORLD.Get_rank()
MPI_SIZE = MPI.COMM_WORLD.Get_size()

gsmh2AlyaCellTypes = {
	1  : 2, #1: 2-node line.
#8: 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
#26: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
#27: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
#28: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
	2  : 10, #2: 3-node triangle.
	3  : 12, #3: 4-node quadrangle.
#9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
#10: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
#16: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
#20: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
#21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
#22: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
#23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
#24: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
#25: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
	4  : 30, #4: 4-node tetrahedron.
	5  : 37, #5: 8-node hexahedron.
	6  : 34, #6: 6-node prism.
	7  : 32, #7: 5-node pyramid.
#11: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
	12 : 39, #12: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
#13: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
#14: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
#17: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
#18: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
#19: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
#29: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
#30: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
#31: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
    36 : 36, #16-node spectral p3 quad CHECK THE ID ALYA, ME LA HE INVENTADO
	92 : 40, #92: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume) 
}

Alya2Names = {
	2  : 'LIN02', #1: 2-node line.
#8: 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
#26: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
#27: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
#28: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
	10 : 'TRI03', #2: 3-node triangle.
	12 : 'QUA04', #3: 4-node quadrangle.
#9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
#10: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
#16: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
#20: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
#21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
#22: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
#23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
#24: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
#25: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
	30 : 'TET04', #4: 4-node tetrahedron.
	37 : 'HEX08', #5: 8-node hexahedron.
	34 : 'PEN06', #6: 6-node prism.
	32 : 'PYR05', #7: 5-node pyramid.
#11: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
	39 : 'HEX27', #12: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
#13: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
#14: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
#17: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
#18: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
#19: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
#29: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
#30: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
#31: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
    36 : 'QUA16', #16-node spectral p3 quad
	40 : 'HEX64', #92: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume) 
}

gsmh2Nodes = {
	1  : 2, #1: 2-node line.
#8: 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
#26: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
#27: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
#28: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
	2  : 3, #2: 3-node triangle.
	3  : 4, #3: 4-node quadrangle.
#9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
#10: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
#16: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
#20: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
#21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
#22: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
#23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
#24: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
#25: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
	4  : 4, #4: 4-node tetrahedron.
	5  : 8, #5: 8-node hexahedron.
	6  : 6, #6: 6-node prism.
	7  : 5, #7: 5-node pyramid.
#11: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
	12 : 27, #12: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
#13: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
#14: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
#17: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
#18: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
#19: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
#29: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
#30: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
#31: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
    36 : 16, #16-node spectral p3 quad
	92 : 64, #92: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume) 
}

def _reduce(v1,v2,dtype):
	'''
	v1 and v2 are two arrays that contain
	v1[:,0] -> index of the minimum
	v1[:,1] -> value of the minimum
	They have both the same size
	'''
	valmat = np.vstack([v1,v2])
	imax   = np.argmax(valmat,axis=0)
	return np.array([valmat[imax[i],i] for i in range(imax.shape[0])],v1.dtype)

gmsh_reduce = MPI.Op.Create(_reduce, commute=True)

def find_node_in_elems(inode,conec):
	'''
	Finds if a node ID is inside a connectivity list
	and returns which element has the node.
	'''
	return np.where(np.any(np.isin(conec,inode),axis=1))[0]

#------------------------------------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------------------------------------

def raiseError(errmsg,all=True):
	'''
	Raise a controlled error and abort execution on
	all processes.
	'''
	if all:
		print('Error: %d - %s' % (MPI_RANK,errmsg),file=sys.stderr,flush=True)
	else:
		if MPI_RANK == 0 or MPI_SIZE == 1:
			print('Error: %d - %s' % (MPI_RANK,errmsg),file=sys.stderr,flush=True)
	MPI_COMM.Abort(1)

def h5_save_field(fname,xyz,varDict,mpio=True,write_master=False,metadata={}):
	'''
	Save a field in HDF5
	'''
	if mpio and not MPI_SIZE == 1:
		h5_save_field_mpio(fname,xyz,varDict,write_master,metadata)
	else:
		h5_save_field_serial(fname,xyz,varDict,metadata)

def h5_save_field_serial(fname,xyz,varDict,metadata={}):
	'''
	Save a field in HDF5 in serial mode
	'''
	# Open file for writing
	file = h5py.File(fname,'w')
	# Metadata
	meta_group = file.create_group('metadata')
	# Store number of points
	dset = meta_group.create_dataset('npoints',(1,),dtype='i')
	dset[:] = xyz.shape[0]
	# Store metadata
	for var in metadata.keys():
		if var in meta_group:
			dset = meta_group[var]
		else:
			dset = meta_group.create_dataset(var,(1,),dtype=metadata[var][1])
		dset[:] = metadata[var][0]
	# Store xyz coordinates
	file.create_dataset('xyz',xyz.shape,dtype=xyz.dtype,data=xyz)
	# Store the variables
	for var in varDict.keys():
		v = varDict[var]
		file.create_dataset(var,v.shape,dtype=v.dtype,data=v)
	file.close()


#------------------------------------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------------------------------------


# Argparse
#pyAlya.cr_start('gmsh2alya',0)
argpar = argparse.ArgumentParser(prog='pyalya_gmsh2alya', description='Export a mesh from GMSH to Alya format')
argpar.add_argument('mesh_file',type=str,help='GMSH input file')
argpar.add_argument('-o','--output',type=str,help='output file name')
argpar.add_argument('-s','--size',type=int,help='maximum size of the data block to read')
argpar.add_argument('-p','--periodic',type=str,help='codes of periodic boundaries, if any, as a string')
argpar.add_argument('--scale',type=str,help='scaling vector (default: 1,1,1)')
argpar.add_argument('-2','--is2D',action='store_true',help='parse a 2D mesh instead')

# Parse inputs
args = argpar.parse_args()
if not args.output:   args.output = args.mesh_file
if not args.periodic: args.periodic = []
else:                 args.periodic = [int(i) for i in args.periodic.split(',')]
if not args.scale:    args.scale    = '1,1,1'

args.scale   = [float(i) for i in args.scale.split(',')]
default_size = True if not args.size else False
dim_id       = 2 if args.is2D else 3



# Print info in screen
print('--|')
print('--| pyalya_gmsh2sod2d |-- ')
print('--|')
print('--| Export a mesh from GMSH to sod2d format.')
print('--|',flush=True)

# Open HDF5 file
h5filename = args.output+'.h5'
h5file = h5py.File(h5filename,'w')
dims_group = h5file.create_group('dims')

# Open GMSH file
args.mesh_file += '.msh'
print('--|')
print('--| Opening <%s>... '%args.mesh_file,end='',flush=True)
mshFile = open(args.mesh_file,'r') #if pyAlya.utils.is_rank_or_serial() else None
print('done!')
print('--|',flush=True)

# Check GMSH version
vers = np.genfromtxt(mshFile,comments='$',max_rows=1) #if #pyAlya.utils.is_rank_or_serial() else None
#vers = #pyAlya.utils.mpi_bcast(vers)
print('--|')
print('--| GMSH file version <%.1f>.'%vers[0])
print('--|',flush=True)
if not vers[0] == 2.2:
	raiseError('This parser can only understand version 2.2 of the Gmsh file format')
# At this point we have checked that the file version is 2.2

# Read the number of zones
nzones = int(np.genfromtxt(mshFile,comments='$',max_rows=1)) #if pyAlya.utils.is_rank_or_serial() else None
#nzones = #pyAlya.utils.mpi_bcast(nzones)
print('--|')
print('--| Detected <%d> physical names. Reading... '%nzones,end='',flush=True)

#if pyAlya.utils.is_rank_or_serial():
# Read from file
data = np.genfromtxt(mshFile,dtype=('i8','i8','<U256'),comments='$',max_rows=nzones)
# Generate a dictionary containing the boundary information
zones = {
	'name'  : np.array([z['f2'].replace('"','') for z in data] if data.ndim > 0 else [data['f2'].tolist().replace('"','')]),
	'code'  : np.array([z['f1'] for z in data] if data.ndim > 0 else [data['f1'].tolist()]),
	'dim'   : np.array([z['f0'] for z in data] if data.ndim > 0 else [data['f0'].tolist()]),
	'isbc'  : np.array([z['f0'] != dim_id for z in data] if data.ndim > 0 else [data['f0'].tolist() != dim_id]),
	'isper' : np.zeros((nzones,),dtype=bool),
}
del data
# Build periodicity
zones['isper'] = [True if zones['code'][iz] in args.periodic else False for iz in range(nzones)]
#else:
#	zones = None
#zones = #pyAlya.utils.mpi_bcast(zones)
print('done!')
print('--|',flush=True)
# Failsafes
if np.sum(np.logical_not(zones['isbc'])) > 1: 
    print('More than one interior zone is not supported!')
	#pyAlya.raiseError('More than one interior zone is not supported!')

# Now read the number of nodes
#pyAlya.cr_start('gmsh2alya.nodes',0)
nnodes = int(np.genfromtxt(mshFile,comments='$',max_rows=1)) #if pyAlya.utils.is_rank_or_serial() else None
#nnodes = pyAlya.utils.mpi_bcast(nnodes)
if default_size: args.size = nnodes
print('--|')
print('--| Detected <%d> nodes.'%nnodes,flush=True)

dset = dims_group.create_dataset('numNodes',(1,),dtype='i8',data=nnodes)

# Generate header for the COORD file
"""
fname  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'COORD')
header = #pyAlya.io.AlyaMPIO_header(
	fieldname   = 'COORD',
	dimension   = 'VECTO',
	association = 'NPOIN',
	npoints     = nnodes,
	nsub        = 1,
	sequence    = 'SEQUE',
	ndims       = dim_id,
	itime       = 0,
	time        = 0,
	ignore_err  = True
)
"""
# Read the number of nodes in batches and write the COORD file
numBatches = int(np.ceil(nnodes/args.size))
print('--| Reading Nodes coords in %d batches of %d...'%(numBatches,args.size),flush=True)

for ibatch in range(numBatches):
	print('--|   Batch %d... '%(ibatch+1),end='',flush=True)
    # Read from text file
	nread = min(args.size,nnodes-ibatch*args.size)
	#print('nread <%d>'%nread,flush=True)
	nodes_data = np.genfromtxt(mshFile,comments='$',max_rows=nread)[:,1:dim_id+1]
	#data  = np.genfromtxt(file,comments='$',max_rows=nread)[:,1:dim_id+1] #if pyAlya.utils.is_rank_or_serial() else None
	# Scale
	for idim in range(dim_id):
		nodes_data[:,idim] *= args.scale[idim]
	if ibatch == 0:
		nodes_dset = h5file.create_dataset('coords',(nread,dim_id),dtype='f8',data=nodes_data,chunks=True,maxshape=(nnodes,dim_id))
	else:
		h5file['coords'].resize((h5file['coords'].shape[0] + nodes_data.shape[0]), axis=0)
		h5file['coords'][-nodes_data.shape[0]:] = nodes_data

	print('done!')
 
print('--|',flush=True)

#nodes_xyz = h5file.create_group('nodes_xyz')
#nodes_dset = h5file.create_dataset('coords',(nread,dim_id),dtype='f8',data=nodes_data)

del nodes_data
#del fname, header
#pyAlya.cr_stop('gmsh2alya.nodes',0)

# Now read the number of elements
#pyAlya.cr_start('gmsh2alya.elements',0)
nelems = int(np.genfromtxt(mshFile,comments='$',max_rows=1)) #if pyAlya.utils.is_rank_or_serial() else None
#nelems = pyAlya.utils.mpi_bcast(nelems)
if default_size: args.size = nelems
print('--|')
print('--| Detected <%d> elements in total.'%nelems,flush=True)
# Boundary and interior elements are stored now consecutively, however,
# at this point we do not know how many of them are present.
#
# We will allocate a memory space for each one of them and we will store
# them as we read by chunks.
#if #pyAlya.utils.is_rank_or_serial():
print('--| Scan for boundary and interior elements.',flush=True)
id_interior  = int(zones['code'][np.logical_not(zones['isbc'])])
nel_interior, nel_boundary = 0, 0
lnods_ndim, lnodb_ndim     = 0, 0
nbatchi, nbatchb           = 0, 0
nel_periodic 					= 0
lnodp_ndim						= 0
nbatchp							= 0

# Read the number of elements in batches 
numBatches = int(np.ceil(nelems/args.size))
print('--| Reading elements in %d batches of %d...'%(numBatches,args.size),flush=True)
for ibatch in range(numBatches):
	print('--|   Batch %d... '%(ibatch+1),end='',flush=True)	
	# Read from text file
	nread = min(args.size,nelems-ibatch*args.size)
	iel_interior, iel_boundary, iel_periodic = 0, 0, 0
	eltyi = -np.ones((nread,),np.int32)
	eltyb = -np.ones((nread,),np.int32)
	eltyp = -np.ones((nread,),np.int32)
	codeb = -np.ones((nread,),np.int32)
	lnods = np.zeros((nread,64),np.int32) # 64 as the maximum order element that we can have
	lnodb = np.zeros((nread,16),np.int32) # 16 as the maximum order element that we can have
	lnodp = np.zeros((nread,16),np.int32) # 16 as the maximum order element that we can have
	# Read the file line by line
	for iline in range(nread):
		# Read one line
		linestr = mshFile.readline()
		line    = np.array([int(l) for l in linestr.split()])
		# Parse line
		eltype = line[1]
		elkind = line[3]
		conec  = line[5:]
		# Skip the element in case of periodicity
		if elkind in args.periodic:
			# Periodic element
			eltyp[iel_periodic]             = gsmh2AlyaCellTypes[eltype]
			lnodp[iel_periodic,:len(conec)] = conec
			#codeb[iel_periodic]             = elkind
			lnodp_ndim                      = max(lnodp_ndim,gsmh2Nodes[eltype])
			iel_periodic 						 += 1
		# We need to find we are interior or boundary
		elif elkind == id_interior:
			# Interior element
			eltyi[iel_interior]             = gsmh2AlyaCellTypes[eltype]
			lnods[iel_interior,:len(conec)] = conec
			lnods_ndim                      = max(lnods_ndim,gsmh2Nodes[eltype])
			iel_interior                   += 1
		else:
			# Boundary element
			eltyb[iel_boundary]             = gsmh2AlyaCellTypes[eltype]
			lnodb[iel_boundary,:len(conec)] = conec
			codeb[iel_boundary]             = elkind
			lnodb_ndim                      = max(lnodb_ndim,gsmh2Nodes[eltype])
			iel_boundary                   += 1
	# Finish the batch read
	nel_interior += iel_interior
	nel_boundary += iel_boundary
	nel_periodic += iel_periodic
	# Get rid of unwanted interior points
	to_keep = eltyi != -1
	eltyi   = eltyi[to_keep]
	lnods   = lnods[to_keep,:]
	# Get rid of unwanted boundary points
	to_keep = eltyb != -1
	eltyb   = eltyb[to_keep]
	lnodb   = lnodb[to_keep,:]
	codeb   = codeb[to_keep]
	# Get rid of unwanted periodic points
	to_keep = eltyp != -1
	eltyp   = eltyp[to_keep]
	lnodp   = lnodp[to_keep,:]
	#codep   = codep[to_keep]
	lnods_len=len(lnods)
	lnodb_len=len(lnodb)
	lnodp_len=len(lnodp)
	lcode_len=len(codeb)
	print('--| lnods_len <%d> lnodb_len <%d> londp_len <%d> lcode_len<%d>'%(lnods_len,lnodb_len,lnodp_len,lcode_len),flush=True)
	if ibatch == 0:
		connec_dset = h5file.create_dataset('connec',(nel_interior,lnods_ndim),dtype='i8',data=lnods,chunks=True,maxshape=(None,lnods_ndim))
		bounds_dset = h5file.create_dataset('boundFaces',(nel_boundary,lnodb_ndim),dtype='i8',data=lnodb,chunks=True,maxshape=(None,lnodb_ndim))
		per_dset = h5file.create_dataset('periodicFaces',(nel_periodic,lnodp_ndim),dtype='i8',data=lnodp,chunks=True,maxshape=(None,lnodp_ndim))
		boundId_dset = h5file.create_dataset('boundFacesId',(nel_boundary,),dtype='i8',data=codeb,chunks=True,maxshape=(None))
	else:
		if lnods_len != 0:
			h5file['connec'].resize((h5file['connec'].shape[0] + lnods.shape[0]), axis=0)
			h5file['connec'][-lnods.shape[0]:] = lnods
		if lnodb_len != 0:
			h5file['boundFaces'].resize((h5file['boundFaces'].shape[0] + lnodb.shape[0]), axis=0)
			h5file['boundFaces'][-lnodb.shape[0]:] = lnodb
		if lnodp_len != 0:
			h5file['periodicFaces'].resize((h5file['periodicFaces'].shape[0] + lnodp.shape[0]), axis=0)
			h5file['periodicFaces'][-lnodp.shape[0]:] = lnodp
		if lcode_len != 0:
			h5file['boundFacesId'].resize((h5file['boundFacesId'].shape[0] + codeb.shape[0]), axis=0)
			h5file['boundFacesId'][-codeb.shape[0]:] = codeb


print('done!')
print('--|',flush=True)
print('--| Found %d inner elements, %d boundary elements, %d per elems'%(nel_interior,nel_boundary,nel_periodic),flush=True)

#lnods_len=len(lnods)
#lnodb_len=len(lnodb)
#lnodp_len=len(lnodp)
#print('--| lnods_len <%d> lnodb_len <%d> londp_len <%d>'%(lnods_len,lnodb_len,lnodp_len),flush=True)

elems_dset = dims_group.create_dataset('numElements',(1,),dtype='i8',data=nel_interior)
bound_dset = dims_group.create_dataset('numBoundaryFaces',(1,),dtype='i8',data=nel_boundary)
per_dset   = dims_group.create_dataset('numPeriodicFaces',(1,),dtype='i8',data=nel_periodic)

#connec_group = h5file.create_group('connectivity')
#connec_dset = h5file.create_dataset('connec',(nel_interior,lnods_ndim),dtype='i4',data=lnods)
#bounds_dset = h5file.create_dataset('boundFaces',(nel_boundary,lnodb_ndim),dtype='i4',data=lnodb)
#per_dset = h5file.create_dataset('periodicFaces',(nel_periodic,lnodp_ndim),dtype='i4',data=lnodp)
#boundId_dset = h5file.create_dataset('boundFacesId',(nel_boundary,),dtype='i4',data=codeb)

del lnods
del lnodb
del lnodp
del codeb

#max_num_per = nel_periodic*16 #de moment hardcodejat 16
dummy_lnodp = np.zeros((0,2),np.int32) 
#---- implemented by me from scratch ----- #
num_bounds_per = 0
if(len(args.periodic) != 0):
	num_bounds_per = np.genfromtxt(mshFile,dtype=('i8'),comments='$',max_rows=1) 
print('--| num bounds periodic: %d'%(num_bounds_per))
npernodes=0
for iper in range(num_bounds_per):
	print('--| Reading per bound %d... '%(iper+1),flush=True)
	linestr = mshFile.readline()
	linestr = mshFile.readline()
	nperlinks = np.genfromtxt(mshFile,dtype=('i8'),comments='$',max_rows=1)
	if default_size: args.size = nperlinks
	print('--| Per bound <%d> per links: %d'%(iper+1,nperlinks),flush=True)
	# Read the number of Per Bounds in batches 
	numBatches = int(np.ceil(nperlinks/args.size))
	print('--| Reading Per Bounds in %d batches of %d...'%(numBatches,args.size),flush=True)
	for ibatch in range(numBatches):
		print('--|   Batch %d... '%(ibatch+1),end='',flush=True)	
		# Read from text file
		nread = min(args.size,nperlinks-ibatch*args.size)

		#data_per = np.genfromtxt(mshFile,comments='$',max_rows=nread)[:,1:dim_id+1]
		data_per = np.genfromtxt(mshFile,comments='$',dtype=('i8'),max_rows=nread)
		if ibatch == 0 and iper == 0:
			nes_dset = h5file.create_dataset('periodicLinks',(nread,2),dtype='i4',data=data_per,chunks=True,maxshape=(None,2))
		else:
			h5file['periodicLinks'].resize((h5file['periodicLinks'].shape[0] + data_per.shape[0]), axis=0)
			h5file['periodicLinks'][-data_per.shape[0]:] = data_per
		del data_per
		print('done!')
	
	npernodes += nperlinks

print('--|',flush=True)
#---- en section implemented by me ------- #

dset = dims_group.create_dataset('numPeriodicLinks',(1,),dtype='i8',data=npernodes)

h5file.close()
mshFile.close()

"""
else:
	nel_interior, nel_boundary = 0, 0
	lnods_ndim, lnodb_ndim     = 0, 0
	nbatchi, nbatchb           = 0, 0
	eltyi, lnods               = None, None
	eltyb, lnodb, codeb        = None, None, None
nel_interior = #pyAlya.utils.mpi_bcast(nel_interior)
nel_boundary = #pyAlya.utils.mpi_bcast(nel_boundary)
lnods_ndim   = #pyAlya.utils.mpi_bcast(lnods_ndim)
lnodb_ndim   = #pyAlya.utils.mpi_bcast(lnodb_ndim)
nbatchi      = #pyAlya.utils.mpi_bcast(nbatchi)
nbatchb      = #pyAlya.utils.mpi_bcast(nbatchb)
#pyAlya.cr_stop('gmsh2alya.elements',0)
"""
# Boundary elements
#pyAlya.cr_start('gmsh2alya.boundary',0)
"""
if nel_boundary > 0:
	print('--| Writing boundary elements.',flush=True)
	# Boundary type
	fname_ltypb  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'LTYPB')
	header_ltypb = #pyAlya.io.AlyaMPIO_header(
		fieldname   = 'LTYPB',
		dimension   = 'SCALA',
		association = 'NBOUN',
		dtype       = 'INTEG',
		size        = '4BYTE',
		npoints     = nel_boundary,
		nsub        = 1,
		sequence    = 'SEQUE',
		ndims       = 1,
		itime       = 0,
		time        = 0,
		ignore_err  = True
	)
	# Boundary connectivity
	fname_lnodb  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'LNODB')
	header_lnodb = #pyAlya.io.AlyaMPIO_header(
		fieldname   = 'LNODB',
		dimension   = 'VECTO',
		association = 'NBOUN',
		dtype       = 'INTEG',
		size        = '4BYTE',
		npoints     = nel_boundary,
		nsub        = 1,
		sequence    = 'SEQUE',
		ndims       = lnodb_ndim,
		itime       = 0,
		time        = 0,
		ignore_err  = True
	)
	# Boundary codes
	fname_codbo  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'CODBO')
	header_codbo = #pyAlya.io.AlyaMPIO_header(
		fieldname   = 'CODBO',
		dimension   = 'SCALA',
		association = 'NBOUN',
		dtype       = 'INTEG',
		size        = '4BYTE',
		npoints     = nel_boundary,
		nsub        = 1,
		sequence    = 'SEQUE',
		ndims       = 1,
		itime       = 0,
		time        = 0,
		ignore_err  = True
	)
	# Boundary sets
	fname_lbset  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'LBSET')
	header_lbset = #pyAlya.io.AlyaMPIO_header(
		fieldname   = 'LBSET',
		dimension   = 'SCALA',
		association = 'NBOUN',
		dtype       = 'INTEG',
		size        = '4BYTE',
		npoints     = nel_boundary,
		nsub        = 1,
		sequence    = 'SEQUE',
		ndims       = 1,
		itime       = 0,
		time        = 0,
		ignore_err  = True
	)

	if default_size:
		# Directly write if the arrays are available on memory
		if #pyAlya.utils.is_rank_or_serial(): lnodb = lnodb[:,:lnodb_ndim]
		# Store into MPIO file
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_ltypb,eltyb,header_ltypb,nel_boundary,0,rank=0)
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lnodb,lnodb,header_lnodb,nel_boundary,0,rank=0)
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_codbo,codeb,header_codbo,nel_boundary,0,rank=0)
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lbset,codeb,header_lbset,nel_boundary,0,rank=0)
	else:
		# We have stored partial files on disk to save memory
		# perform a read by batches and store into MPIO
		nwrite, nskip = 0, 0
		#pyAlya.pprint(0,'--| Writing in batches: ',flush=True)
		for ibatch in range(nbatchb):
			#pyAlya.pprint(0,'--|   Batch %d... '%(ibatch+1),end='',flush=True)
			# Load batch data
			if #pyAlya.utils.is_rank_or_serial():
				data   = np.load('boundary_%d.npz'%(ibatch+1))
				eltyb  = data['eltyb']
				lnodb  = data['lnodb'][:,:lnodb_ndim]
				codeb  = data['codeb']
				nwrite = eltyb.shape[0]
				os.remove('boundary_%d.npz'%(ibatch+1))
			else:
				eltyb = None
				lnodb = None
				codeb = None
			# Store into MPIO file
			#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_ltypb,eltyb,header_ltypb,nwrite,nskip,rank=0)
			#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lnodb,lnodb,header_lnodb,nwrite,nskip,rank=0)
			#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_codbo,codeb,header_codbo,nwrite,nskip,rank=0)
			#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lbset,codeb,header_lbset,nwrite,nskip,rank=0)
			if #pyAlya.utils.is_rank_or_serial(): nskip += nwrite
			#pyAlya.pprint(0,'done!')
	#pyAlya.pprint(0,'--|',flush=True)
	del fname_ltypb, fname_codbo, fname_lbset
	del header_ltypb, header_lnodb, header_codbo, header_lbset
#pyAlya.cr_stop('gmsh2alya.boundary',0)

# Interior elements
#pyAlya.cr_start('gmsh2alya.interior',0)
#pyAlya.pprint(0,'--| Writing interior elements.',flush=True)
# Element type
fname_ltype  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'LTYPE')
header_ltype = #pyAlya.io.AlyaMPIO_header(
	fieldname   = 'LTYPE',
	dimension   = 'SCALA',
	association = 'NELEM',
	dtype       = 'INTEG',
	size        = '4BYTE',
	npoints     = nel_interior,
	nsub        = 1,
	sequence    = 'SEQUE',
	ndims       = 1,
	itime       = 0,
	time        = 0,
	ignore_err  = True
)
# Element connectivity
fname_lnods  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'LNODS')
header_lnods = #pyAlya.io.AlyaMPIO_header(
	fieldname   = 'LNODS',
	dimension   = 'VECTO',
	association = 'NELEM',
	dtype       = 'INTEG',
	size        = '4BYTE',
	npoints     = nel_interior,
	nsub        = 1,
	sequence    = 'SEQUE',
	ndims       = lnods_ndim,
	itime       = 0,
	time        = 0,
	ignore_err  = True
)
eltype = []
if default_size:
	# Directly write if the arrays are available on memory
	if #pyAlya.utils.is_rank_or_serial(): lnods = lnods[:,:lnods_ndim]
	# Store into MPIO file
	#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_ltype,eltyi,header_ltype,nel_interior,0,rank=0)
	#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lnods,lnods,header_lnods,nel_interior,0,rank=0)
	if #pyAlya.utils.is_rank_or_serial(): eltype = [Alya2Names[e] for e in np.unique(eltyi)]
else:
	# We have stored partial files on disk to save memory
	# perform a read by batches and store into MPIO
	nwrite, nskip = 0, 0
	#pyAlya.pprint(0,'--| Writing in batches: ',flush=True)
	for ibatch in range(nbatchi):
		#pyAlya.pprint(0,'--|   Batch %d... '%(ibatch+1),end='',flush=True)
		# Load batch data
		if #pyAlya.utils.is_rank_or_serial():
			data   = np.load('interior_%d.npz'%(ibatch+1))
			eltyi  = data['eltyi']
			lnods  = data['lnods'][:,:lnods_ndim]
			nwrite = eltyi.shape[0]
			os.remove('interior_%d.npz'%(ibatch+1))
			eltype = np.hstack([eltype,np.unique(eltyi)])
		else:
			eltyi = None
			lnods = None
		# Store into MPIO file
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_ltype,eltyi,header_ltype,nwrite,nskip,rank=0)
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lnods,lnods,header_lnods,nwrite,nskip,rank=0)
		if #pyAlya.utils.is_rank_or_serial(): nskip += nwrite
		#pyAlya.pprint(0,'done!')
	if #pyAlya.utils.is_rank_or_serial(): eltype = [Alya2Names[e] for e in np.unique(eltyi)]
#pyAlya.pprint(0,'--|',flush=True)
del fname_ltype, header_ltype, header_lnods
#pyAlya.cr_stop('gmsh2alya.interior',0)

# Compute LELBO
# Idea 2, read by zones both connectivities and find the faces
# This is a real pain in the ass but it must be done
#pyAlya.utils.mpi_barrier()
#pyAlya.cr_start('gmsh2alya.lelbo',0)
if nel_boundary > 0:
	#pyAlya.pprint(0,'--| Computing element to boundary relation.')
	fname_lelbo  = #pyAlya.io.MPIO_AUXFILE_S_FMT % (args.casename,'LELBO')
	header_lelbo = #pyAlya.io.AlyaMPIO_header(
		fieldname   = 'LELBO',
		dimension   = 'SCALA',
		association = 'NBOUN',
		dtype       = 'INTEG',
		size        = '4BYTE',
		npoints     = nel_boundary,
		nsub        = 1,
		sequence    = 'SEQUE',
		ndims       = 1,
		itime       = 0,
		time        = 0,
		ignore_err  = True
	)
	if default_size:
		# Assume LNODS to be already in memory
		lnods = #pyAlya.utils.mpi_scatter(lnods,root=0,do_split=True) # scatter LNODS
	else:
		# Compute a worksplit and have each processor load a part
		# of the element connectivity matrix
		mystart,myend = #pyAlya.utils.worksplit(0,nel_interior,MPI_RANK)
		lnods,_       = #pyAlya.io.AlyaMPIO_readByChunk(fname_lnods,myend-mystart,mystart)
	# Obtain the lenghts of lnods for the offset computation
	elemlist = np.arange(lnods.shape[0],dtype=np.int32)
	lens     = #pyAlya.utils.mpi_gather(lnods.shape[0],all=True)
	offst    = max(np.sum(lens[:MPI_RANK]),0) if MPI_SIZE > 1 else 0

	# Now compute and write LELBO
	if default_size: args.size = nel_boundary
	#pyAlya.pprint(0,'--| Writing element to boundary relation in batches of %d: '%args.size,flush=True)
	for ibatch in range(int(np.ceil(nel_boundary/args.size))):
		#pyAlya.pprint(0,'--|   Batch %d... '%(ibatch+1),end='',flush=True)
		rows2read = min(args.size,nel_boundary-ibatch*args.size)
		rows2skip = ibatch*args.size
		# All processors should read the batch
		if default_size:
			lnodb   = #pyAlya.utils.mpi_bcast(lnodb,root=0) # scatter LNODB
		else:
			lnodb,_ = #pyAlya.io.AlyaMPIO_readByChunk(fname_lnodb,rows2read,rows2skip)
		vals = -np.ones((rows2read,),np.int32)
		# Loop all the elements in the batch
		for ib in range(rows2read):
			# Find all the interior elements that contain the first node
			ielems = find_node_in_elems(lnodb[ib,0],lnods)
			candidate_elements = elemlist[ielems]
			# Now, we have one of the following 3 options:
			#   1. We have multiple candidates, then we need to keep scanning nodes
			#   2. We have one candidate, then we have found the element
			#   3. We have found no candidate, hence we don't have the element
			lnods2    = lnods[ielems,:]
			elemlist2 = elemlist[ielems]
			# Loop the remaining nodes and find candidates
			for ii in range(1,lnodb.shape[1]):
				if lnodb[ib,ii] == 0: continue
				# Search on a reduced list
				ielems    = find_node_in_elems(lnodb[ib,ii],lnods2)
				candidate_elements = elemlist2[ielems]
				lnods2    = lnods2[ielems,:]
				elemlist2 = elemlist2[ielems]
			# Someone will have found the element, write it to vals
			if len(candidate_elements) == 1: vals[ib] = candidate_elements[0] + offst + 1
		# Reduce
		vals = #pyAlya.utils.mpi_reduce(vals,op=gmsh_reduce,all=True)
		# Store into MPIO file
		#pyAlya.io.AlyaMPIO_writeByChunk_serial(fname_lelbo,vals,header_lelbo,rows2read,rows2skip,rank=0)
		#pyAlya.pprint(0,'done!',flush=True)
	#pyAlya.pprint(0,'--|',flush=True)
	del fname_lelbo, fname_lnods, fname_lnodb, header_lelbo
#pyAlya.cr_stop('gmsh2alya.lelbo',0)

# Finishing touches
#pyAlya.pprint(0,'--| Finishing mesh export...',flush=True)

# Create .dims.dat file
#if pyAlya.utils.is_rank_or_serial():
	filedims = open(args.casename+'.dims.dat','w')
	filedims.write('NODAL_POINTS            %d \n'%nnodes)
	filedims.write('ELEMENTS                %d \n'%nel_interior)
	filedims.write('SPACE_DIMENSIONS        %d \n'%dim_id)
	filedims.write('TYPES_OF_ELEMENTS       %s'%eltype[0])
	for e in eltype[1:]: filedims.write(',%s'%e)
	filedims.write('\n')
	filedims.write('BOUNDARIES              %d \n'%nel_boundary)
	filedims.close()
#pyAlya.pprint(0,'--| Written dimensions file <%s>.'%(args.casename+'.dims.dat'),flush=True)

# Write info file
#if pyAlya.utils.is_rank_or_serial():
	fileinfo = open(args.casename+'.info.dat','w')
	fileinfo.write('BOUNDARY_CODES\n')
	for iz in range(nzones):
		if zones['isbc'][iz] and not zones['isper'][iz]:
			fileinfo.write('    %d    %s\n'%(zones['code'][iz],zones['name'][iz]))
	fileinfo.write('END_BOUNDARY_CODES\n')
	fileinfo.close()
#pyAlya.pprint(0,'--| Written info file <%s>.'%(args.casename+'.dims.dat'),flush=True)
#pyAlya.pprint(0,'--|',flush=True)

# Close and say goodbye
#if pyAlya.utils.is_rank_or_serial(): file.close()
#pyAlya.pprint(0,'--|')
#pyAlya.pprint(0,'--| Bye!',flush=True)
#pyAlya.cr_stop('gmsh2alya',0)
#pyAlya.cr_info()
"""
