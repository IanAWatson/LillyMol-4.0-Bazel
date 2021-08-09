#ifndef MOLECULE_LIB_COORDINATE_BOX_H
#define MOLECULE_LIB_COORDINATE_BOX_H

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/space_vector.h"

namespace coordinate_box {

// The idea is that by establishing a box around a molecule, and a grid,
// then the position of any atom can be specified by a single number. As long
// as there is a means of transforming from that single number to a cell
// within the box.
// Since this is 3 dimensions, there are two dimensions for the box - it
// is unconstrainted in the z dimension.
// Set CELL_NUMBER_CHECK_IN_BOUNDS to enable bounds checking.
class CoordinateBox {
  private:
    // The length of each cell.
    double _cell_size;

    // The number of cells along each side of the box. Z is unspecified.
    int _cell[2];

    // Once a position in the box is determined, that then gets
    // translated to where the origin is.
    Space_Vector<double> _origin;

  public:
    CoordinateBox();

    // Syntax is %B{cell_size,x,y,z}. This function assumes
    // just the csv values.
    int BuildFromSmilesToken(const const_IWSubstring& token);

    // For a position in space, the corresponding cell.
    template <typename T>
    int CellNumber(const Space_Vector<T>& coords) const {
      return CellNumber(coords.x(), coords.y(), coords.z());
    }
    template <typename T>
    int CellNumber(T x, T y, T z) const;


    // For a cell number, return the position in space.
    template <typename T> 
    std::tuple<T, T, T> CoordinatesAsTuple(int cell_number) const;
    template <typename T>
    Space_Vector<T> CoordinatesAsVector(int cell_number) const;
};

} // namespace coordinate_box

#endif // MOLECULE_LIB_COORDINATE_BOX_H
