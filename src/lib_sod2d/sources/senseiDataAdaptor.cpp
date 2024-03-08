#include "senseiDataAdaptor.h"

#include <mpi.h>
#include <svtkCellArray.h>
#include <svtkCellData.h>
#include <svtkDoubleArray.h>
#include <svtkFloatArray.h>
#include <svtkIdTypeArray.h>
#include <svtkImageData.h>
#include <svtkIntArray.h>
#include <svtkNew.h>
#include <svtkObjectFactory.h>
#include <svtkPointData.h>
#include <svtkPoints.h>
#include <svtkPolyData.h>
#include <svtkRectilinearGrid.h>
#include <svtkSmartPointer.h>
#include <svtkUnsignedCharArray.h>
#include <svtkUnstructuredGrid.h>
#include <svtkLagrangeHexahedron.h>
#include <svtkCompositeDataSet.h>
#include <svtkSOADataArrayTemplate.h>
// #include <svtkXMLUnstructuredGridWriter.h>

senseiNewMacro(DataAdaptor);

DataAdaptor::~DataAdaptor()
{
  setGrid(nullptr);
}

int DataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes)
{
  numMeshes = 1;
  return 0;
}

size_t getNumComponents(const std::array<double *, 3> &array)
{
  return 3 - std::count(array.begin(), array.end(), nullptr);
}

int DataAdaptor::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata)
{
  if (id > 0)
  {
    SENSEI_ERROR("invalid mesh id " << id)
    return -1;
  }

  metadata->MeshName = "mesh";

  metadata->MeshType = SVTK_UNSTRUCTURED_GRID;
  metadata->CoordinateType = SVTK_FLOAT;
  metadata->NumBlocks = m_numBlocks;
  metadata->NumGhostCells = 0;
  metadata->NumArrays = 1;
  metadata->StaticMesh = 1;

  for (const auto &field : m_fields)
  {
    metadata->ArrayName.push_back(field.first);
    metadata->ArrayCentering.push_back(svtkDataObject::POINT);
    metadata->ArrayComponents.push_back(3); //xyz
    metadata->ArrayType.push_back(SVTK_FLOAT);
    metadata->NumBlocksLocal.push_back(1);
  }
  return 0;
}

int DataAdaptor::GetMesh(const std::string &meshName, bool structureOnly, svtkDataObject *&mesh)
{

  if(m_grid)
  {
    mesh = m_grid;
    return 0;
  }  
  auto numPointsPerElement = 64; //4x4x4 SVTK_LAGRANGE_HEXAHEDRON
  auto numElements = numConnections / numPointsPerElement;
  auto grid = svtk::MakeSmartPointer(svtkUnstructuredGrid::New());
  std::cerr << "creating mesh with " << numVertices << " verticees, " << numElements << " elements and " << numConnections << " connections" << std::endl;
  auto points = svtk::MakeSmartPointer(svtkPoints::New());
  auto pointsArray = svtk::MakeSmartPointer(svtkSOADataArrayTemplate<float>::New());
  pointsArray->SetNumberOfComponents(3);
  pointsArray->SetArray(0, xyz, numVertices, 1);
  pointsArray->SetArray(1, xyz + numVertices, numVertices, 1);
  pointsArray->SetArray(2, xyz + 2*numVertices, numVertices, 1);
  
  points->SetData(pointsArray);
  grid->SetPoints(points);

  grid->Allocate(numElements);
  for (size_t i = 0; i < numElements; i++)
  {
      
      auto l = svtk::MakeSmartPointer(svtkIdList::New());
      l->SetNumberOfIds(numPointsPerElement);
      for (size_t j = 0; j < numPointsPerElement; j++)
      {
          l->SetId(j, connectivity[i*numPointsPerElement + j]);
      }
      grid->InsertNextCell(SVTK_LAGRANGE_HEXAHEDRON, l);
  }
  
  mesh = grid;
  setGrid(grid);

  return 0;
}

int DataAdaptor::AddArray(svtkDataObject *mesh, const std::string &meshName, int association, const std::string &arrayName)
{
  auto array = m_fields.find(arrayName);
  if (array != m_fields.end())
  {
    svtkUnstructuredGrid *ugrid = dynamic_cast<svtkUnstructuredGrid *>(mesh);
    if (!ugrid)
    {
      if(m_mpiRank == 0)
        std::cerr << meshName << " is not an unstructured grid, can only add arrays to unstructured grids, mesh " << std::endl; 
      return -1;
    }
    auto field = svtk::MakeSmartPointer(svtkSOADataArrayTemplate<float>::New());
    field->SetNumberOfComponents(3);
    
    field->SetArray(0, array->second, numVertices, 1);
    field->SetArray(1, array->second + numVertices, numVertices, 1);
    field->SetArray(2, array->second + 2*numVertices, numVertices, 1);

    field->SetName(arrayName.c_str());
    ugrid->GetPointData()->AddArray(field);
    m_fields.erase(array);
    if(m_mpiRank == 0)
      std::cerr << "added array " << arrayName << " to mesh " << meshName << std::endl;
    return 0;
  }
  if(m_mpiRank == 0)
    std::cerr << "no array named " << arrayName << " on mesh " << meshName << std::endl;
  return -1;
}

int DataAdaptor::AddGhostCellsArray(svtkDataObject *mesh, const std::string &meshName)
{
  return 0;
}

int DataAdaptor::ReleaseData()
{
  std::cerr << "ReleaseData called" << std::endl;
  return 0;
}

void DataAdaptor::setGrid(svtkUnstructuredGrid *grid)
{
  if (m_grid && grid != m_grid)
    m_grid->Delete();
  m_grid = grid;
  m_createdGrid = true;
}

void DataAdaptor::setGrid(float *xyz, const int *numVertices, int* connectivity, const int *numConnections)
{
  this->xyz = xyz;
  this->numVertices = *numVertices;
  this->connectivity = connectivity;
  this->numConnections = *numConnections;
}

// void DataAdaptor::add_scalar_field_(double *data, char *name)
// {
//   m_fields[name] = {data, nullptr, nullptr};
// }

void DataAdaptor::add_vector_field_(float *xyzdata, char *name)
{
  m_fields[name] = xyzdata;
}

void DataAdaptor::setNumBlocks(int num)
{
  m_numBlocks = num;
}

void DataAdaptor::setMpiRank(int rank)
{
  m_mpiRank = rank;
}

