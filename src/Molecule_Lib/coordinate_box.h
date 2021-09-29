#ifndef MOLECULE_LIB_COORDINATE_BOX_H
#define MOLECULE_LIB_COORDINATE_BOX_H

#include <cstdint>

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

// A cell number that defines a set of coordinates in relation to a
// concentric set of cubic shells. The shells are defined by a resolution
// specified in the constructor.
// The cell numbering is not public.
// Points going through the conversion process will be truncated at
// values that are multiples of `resolution`.
class ConcentricBox {
  private:
    double _resolution;
    // Several computations need _resolution*0.5
    double _half_resolution;

  // private functions.
    template <typename T> Space_Vector<T> _right_face(uint32_t layer, int within_layer) const;
    template <typename T> Space_Vector<T> _left_face(uint32_t layer, int within_layer) const;
    template <typename T> Space_Vector<T> _somewhere_along_x(uint32_t layer, int within_layer) const;

    template <typename T> uint32_t CellNumber(const Space_Vector<T>& coords, uint32_t layer, uint32_t x, uint32_t y, uint32_t z) const;
    template <typename T> uint32_t LayerNumber(const Space_Vector<T>& coords) const;

  public:
    // A default resolution which should be good enough for most interatomic distances.
    ConcentricBox(double resolution = 0.001);

    // Given a position in space, return the corresponding call number.
    template <typename T>
    uint32_t CellNumber(const Space_Vector<T>& coords) const;

    // For a cell number, return the position in space.
    template <typename T> 
    std::tuple<T, T, T> CoordinatesAsTuple(uint32_t cell_number) const;
    template <typename T>
    Space_Vector<T> CellNumberToCoordinates(uint32_t cell_number) const;
};

// Helper function for ConcentricBox, exposed just for testing.
uint32_t CellToLayer(uint32_t cell);

} // namespace coordinate_box

#endif // MOLECULE_LIB_COORDINATE_BOX_H
