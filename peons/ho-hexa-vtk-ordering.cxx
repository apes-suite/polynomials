// Code to find the ordering index given an ijk tuple and the
// polynomial order for High Order Hexahedron.
//
// source: https://gitlab.kitware.com/vtk/vtk/-/blob/master/Common/DataModel/vtkHigherOrderHexahedron.cxx#L601
//
// Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// BSD-3-Clause
//
// This file can be compiled into an executable that will print
// a sample high-order hexahedron ordering as a Fortran module
// to standard out and can be used to generate a Fortran source
// file providing that sample for testing:
// g++ ho-hexa-vtk-ordering.cxx
// peons/a.out > utests/sample_vtkHOhex_module.f90

#include <iostream>

int PointIndexFromIJK(int i, int j, int k, const int* order)
{
  bool ibdy = (i == 0 || i == order[0]);
  bool jbdy = (j == 0 || j == order[1]);
  bool kbdy = (k == 0 || k == order[2]);
  // How many boundaries do we lie on at once?
  int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

  if (nbdy == 3) // Vertex DOF
  {              // ijk is a corner node. Return the proper index (somewhere in [0,7]):
    return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
  }

  int offset = 8;
  if (nbdy == 2) // Edge DOF
  {
    if (!ibdy)
    { // On i axis
      return (i - 1) + (j ? order[0] + order[1] - 2 : 0) + (k ? 2 * (order[0] + order[1] - 2) : 0) +
        offset;
    }
    if (!jbdy)
    { // On j axis
      return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
        (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
    }
    // !kbdy, On k axis
    offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
    return (k - 1) + (order[2] - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
  }

  offset += 4 * (order[0] + order[1] + order[2] - 3);
  if (nbdy == 1) // Face DOF
  {
    if (ibdy) // On i-normal face
    {
      return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) +
        offset;
    }
    offset += 2 * (order[1] - 1) * (order[2] - 1);
    if (jbdy) // On j-normal face
    {
      return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) +
        offset;
    }
    offset += 2 * (order[2] - 1) * (order[0] - 1);
    // kbdy, On k-normal face
    return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) +
      offset;
  }

  // nbdy == 0: Body DOF
  offset += 2 *
    ((order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) +
      (order[0] - 1) * (order[1] - 1));
  return offset + (i - 1) + (order[0] - 1) * ((j - 1) + (order[1] - 1) * ((k - 1)));
}

// Main program to print the index mapping of an ijk ordering
// and the VTK high-order hexahedron ordering.
//int main(int argc, char* argv[]) {
int main() {
    int order[3] = {9, 8, 7};
    std::cout << "! Created by peons/ho-hexa-vtk-ordering.cxx\n";
    std::cout << "module sample_vtkHOhex_module\n";
    std::cout << "  implicit none\n";
    std::cout << "  private\n";
    std::cout << "  integer, public :: hexorders(3) = [" << order[0]
                                                      << ", " << order[1]
                                                      << ", " << order[2] << "]\n";
    std::cout << "  integer, public :: indexMap(0:" << order[0] << ", 0:" << order[1]
                                                    << ", 0:" << order[2] << ")\n";
    std::cout << "  public :: fillMap\n";
    std::cout << "\n\n";
    std::cout << "contains\n";
    std::cout << "\n\n";

    std::cout << "  subroutine fillMap()\n";
    for (int k=0; k <= order[2]; ++k) {
        for (int j=0; j <= order[1]; ++j) {
            for (int i=0; i <= order[0]; ++i) {
                std::cout << "    indexMap(" << i << ", " << j << ", " << k << ") = "
                          << PointIndexFromIJK(i, j, k, order) << "\n";
            }
        }
    }
    std::cout << "  end subroutine fillMap\n";
    std::cout << "end module sample_vtkHOhex_module\n";
    return 0;
}
