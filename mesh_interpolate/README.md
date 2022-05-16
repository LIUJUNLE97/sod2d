# About

Simple tool to interpolate results from the linear mesh into a high-order one. Assumes the high-order mesh is obtained from the linear one, such that element corner nodes match in both coordinates and index.

# Usage

Generate a linear mesh with GMSH, then populate it with the high-order nodes, storing both versions. The program will read the field values obtained on the linear mesh and use the high-order node locations to interpolate, generating the corresponding high-order field. It then outputs into a sod2d readable format.
