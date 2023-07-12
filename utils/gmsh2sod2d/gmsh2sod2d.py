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

MPI_RANK = MPI.COMM_WORLD.Get_rank()
MPI_SIZE = MPI.COMM_WORLD.Get_size()


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
argpar = argparse.ArgumentParser(prog='gmsh2sod2d script', description='Export a mesh from GMSH to sod2d mesh conversor tool input format')
argpar.add_argument('mesh_file',type=str,help='GMSH input file')
argpar.add_argument('-o','--output',type=str,help='output file name')
argpar.add_argument('-r','--order',type=int,help='mesh order')
argpar.add_argument('-s','--size',type=int,help='maximum size of the data block to read')
argpar.add_argument('-p','--periodic',type=str,help='codes of periodic boundaries, if any, as a string')
argpar.add_argument('--scale',type=str,help='scaling vector (default: 1,1,1)')
argpar.add_argument('-2','--is2D',action='store_true',help='parse a 2D mesh instead')

# Parse inputs
args = argpar.parse_args()
if not args.output:   args.output = args.mesh_file
if not args.periodic: args.periodic = []
else:                 args.periodic = [int(i) for i in args.periodic.split(',')]
if not args.order:    args.order 	= 3
if not args.scale:    args.scale    = '1,1,1'

args.scale   = [float(i) for i in args.scale.split(',')]
default_size = True if not args.size else False
dim_id       = 2 if args.is2D else 3



# Print info in screen
print('--|')
print('--| gmsh2sod2d |-- ')
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
mshFile = open(args.mesh_file,'r') 
print('done!')
print('--|',flush=True)

# Check GMSH version
vers = np.genfromtxt(mshFile,comments='$',max_rows=1)
print('--|')
print('--| GMSH file version <%.1f>.'%vers[0])
print('--|',flush=True)
if not vers[0] == 2.2:
	raiseError('This parser can only understand version 2.2 of the Gmsh file format')
# At this point we have checked that the file version is 2.2

# Order of the mesh
porder = args.order
nnode = (porder+1)**3
npbou = (porder+1)**2
print('--|')
print('--| Order: <%d>, nnode: <%d>, npbou: <%d> '%(porder,nnode,npbou),end='',flush=True)

dset = dims_group.create_dataset('order',(1,),dtype='i8',data=args.order)

# Read the number of zones
nzones = int(np.genfromtxt(mshFile,comments='$',max_rows=1)) 
print('--|')
print('--| Detected <%d> physical names. Reading... '%nzones,end='',flush=True)

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

print('done!')
print('--|',flush=True)
# Failsafes
if np.sum(np.logical_not(zones['isbc'])) > 1: 
    print('More than one interior zone is not supported!')

# Now read the number of nodes

nnodes = int(np.genfromtxt(mshFile,comments='$',max_rows=1)) 

if default_size: args.size = nnodes
print('--|')
print('--| Detected <%d> nodes.'%nnodes,flush=True)

dset = dims_group.create_dataset('numNodes',(1,),dtype='i8',data=nnodes)

# Read the number of nodes in batches and write the COORD file
numBatches = int(np.ceil(nnodes/args.size))
print('--| Reading Nodes coords in %d batches of %d...'%(numBatches,args.size),flush=True)

for ibatch in range(numBatches):
	print('--|   Batch %d... '%(ibatch+1),end='',flush=True)
    # Read from text file
	nread = min(args.size,nnodes-ibatch*args.size)
	#print('nread <%d>'%nread,flush=True)
	nodes_data = np.genfromtxt(mshFile,comments='$',max_rows=nread)[:,1:dim_id+1]
	#data  = np.genfromtxt(file,comments='$',max_rows=nread)[:,1:dim_id+1] 
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

del nodes_data

# Now read the number of elements
nelems = int(np.genfromtxt(mshFile,comments='$',max_rows=1)) 

if default_size: args.size = nelems
print('--|')
print('--| Detected <%d> elements in total.'%nelems,flush=True)
# Boundary and interior elements are stored now consecutively, however,
# at this point we do not know how many of them are present.
#
# We will allocate a memory space for each one of them and we will store
# them as we read by chunks.

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
	lnods = np.zeros((nread,nnode),np.int32)
	lnodb = np.zeros((nread,npbou),np.int32)
	lnodp = np.zeros((nread,npbou),np.int32)
	# Read the file line by line
	for iline in range(nread):
		# Read one line
		linestr = mshFile.readline()
		line    = np.array([int(l) for l in linestr.split()])
		# Parse line
		eltype = line[1]
		elkind = line[3]
		conec  = line[5:]
		#print('eltype <%d> len(connec) <%d>'%(eltype,len(conec)),flush=True)
		# Skip the element in case of periodicity
		if elkind in args.periodic:
			# Periodic element
			if len(conec)!=npbou:
				raiseError('Error! len(conec)!=npbou in periodic element!!')
			eltyp[iel_periodic]             = eltype
			lnodp[iel_periodic,:len(conec)] = conec
			lnodp_ndim                      = max(lnodp_ndim,npbou)
			iel_periodic 						 += 1
		# We need to find we are interior or boundary
		elif elkind == id_interior:
			# Interior element
			if len(conec)!=nnode:
				raiseError('Error! len(conec)!=nnode in interior element!!')
			eltyi[iel_interior]             = eltype
			lnods[iel_interior,:len(conec)] = conec
			lnods_ndim                      = max(lnods_ndim,nnode)
			iel_interior                   += 1
		else:
			# Boundary element
			if len(conec)!=npbou:
				raiseError('Error! len(conec)!=npbou in boundary element!!')
			eltyb[iel_boundary]             = eltype
			lnodb[iel_boundary,:len(conec)] = conec
			codeb[iel_boundary]             = elkind
			lnodb_ndim                      = max(lnodb_ndim,npbou)
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
			nes_dset = h5file.create_dataset('periodicLinks',(nread,2),dtype='i8',data=data_per,chunks=True,maxshape=(None,2))
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
