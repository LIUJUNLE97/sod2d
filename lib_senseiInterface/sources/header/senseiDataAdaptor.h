#pragma once

#include <DataAdaptor.h>
#include <string>
#include <map>

class svtkDoubleArray;
class svtkUnstructuredGrid;
class DataAdaptor : public sensei::DataAdaptor
{
public:
  static DataAdaptor *New();
  ~DataAdaptor();
  senseiTypeMacro(DataAdaptor, sensei::DataAdaptor);
  // SENSEI API : return 0 if success
  int GetNumberOfMeshes(unsigned int &numMeshes) override;

  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &md) override;

  int GetMesh(const std::string &meshName, bool structureOnly,
              svtkDataObject *&mesh) override;

  int AddArray(svtkDataObject *mesh, const std::string &meshName,
               int association, const std::string &arrayName) override;

  int AddGhostCellsArray(svtkDataObject *mesh, const std::string &meshName) override;

  int ReleaseData() override;

  void setNumBlocks(int num);
  void setMpiRank(int rank);

  void setGrid(float *xyz, const int *numVertices, int* connectivity, const int *numConnections);
  // void add_scalar_field_(double *data, char *name);
  void add_vector_field_(float *xyzdata, char *name);
private:
  std::map<std::string, float *> m_fields;
  svtkUnstructuredGrid *m_grid = nullptr;
  // grid data
  float *xyz = nullptr;
  int numVertices = 0;
  int *connectivity = nullptr;
  int numConnections = 0;

  int dim = 3;
  bool m_createdGrid = false;
  int m_numBlocks;
  int m_mpiRank = 0;
  void setGrid(svtkUnstructuredGrid *grid);
};