
extern "C" void print_sensei_(const char *text);
extern "C" void print_sensei_int_(const int *text);

extern "C" void testcoproc_(int *flag);

extern "C" void init_sensei_(const int *mpiRank, const int *mpiSize);

extern "C" void finalize_sensei_();

extern "C" void process_sensei_(int *step);

// this only passes the data required to create the grid to the adaptor
// called creategrid to be conform with catalyst
// hexahedron elements 64 nodes / 4 points per edge
extern "C" void creategrid_(float *xyz, const int *numVertices, int* connectivity, const int *numConnections);

// extern "C" void add_scalar_field_(const double *data, char *name);

extern "C" void add_vector_field_(float *xyzdata, char *name);
