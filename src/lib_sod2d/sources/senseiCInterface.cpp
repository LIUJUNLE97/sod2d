#include "senseiCInterface.h"
#include "senseiDataAdaptor.h"
#include <ConfigurableAnalysis.h>
#include <svtkCellArray.h>
#include <svtkDoubleArray.h>
#include <svtkIdTypeArray.h>
#include <svtkPointData.h>
#include <svtkPoints.h>
#include <svtkUnsignedCharArray.h>
#include <svtkUnstructuredGrid.h>

#include <iostream>

svtkSmartPointer<sensei::ConfigurableAnalysis> analysisAdaptor;
svtkSmartPointer<DataAdaptor> dataAdaptor;

extern "C" void testcoproc_(int *flag)
{
  *flag = analysisAdaptor != nullptr && dataAdaptor != nullptr;
}

void print_sensei_(const char *text)
{
  std::cerr << text << std::endl;
}

void print_sensei_int_(const int *text)
{
  std::cerr << *text << std::endl;
}


void init_sensei_(const int *mpiRank, const int *mpiSize)
{
  analysisAdaptor = svtkSmartPointer<sensei::ConfigurableAnalysis>::New();
  std::cerr << "using sensei config " << SENSEI_CONFIG_FILE << std::endl;
  if (analysisAdaptor->Initialize(SENSEI_CONFIG_FILE))
  {
    std::cerr << "Failed to initialize the analysis adaptor" << std::endl;
  }
  dataAdaptor = svtkSmartPointer<DataAdaptor>::New();
  dataAdaptor->setNumBlocks(*mpiSize);
  dataAdaptor->setMpiRank(*mpiRank);
}

void finalize_sensei_()
{
  if (analysisAdaptor)
    analysisAdaptor->Finalize();
  analysisAdaptor = nullptr;
  dataAdaptor = nullptr;
}

void process_sensei_(int *step)
{
  // std::cerr << "processing alya time step " << *step << std::endl;
  if (analysisAdaptor && dataAdaptor)
  {
    dataAdaptor->SetDataTimeStep(*step);
    // dataAdaptor->SetDataTime(*time);
    sensei::DataAdaptor* dummy;
    analysisAdaptor->Execute(dataAdaptor, &dummy);
  }
  else
  {
    std::cerr << "no analysisAdaptor && dataAdaptor" << std::endl;
  }
}

void creategrid_(float *xyz, const int *numVertices, int* connectivity, const int *numConnections)
{
  dataAdaptor->setGrid(xyz, numVertices, connectivity, numConnections);
}

// void add_scalar_field_(double *data, char *name)
// {
//   dataAdaptor->add_scalar_field_(data, name);
// }

void add_vector_field_(float *xyzData, char *name)
{
  dataAdaptor->add_vector_field_(xyzData, name);
}
