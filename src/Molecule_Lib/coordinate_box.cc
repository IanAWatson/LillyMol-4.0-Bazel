#include "coordinate_box.h"

namespace coordinate_box {
using std::cerr;

CoordinateBox::CoordinateBox() {
  _cell_size = 0.0;
  _cell[0] = 0.0;
  _cell[1] = 0.0;
}

int
CoordinateBox::BuildFromSmilesToken(const const_IWSubstring& from_smiles) {
  resizable_array_p<const_IWSubstring> tokens;
  if (from_smiles.split(tokens, ',') != 3) {
    cerr << "CoordinateBox::BuildFromSmilesToken:must contain 3 tokens " << from_smiles << '\n';
    return 0;
  }

  for (const const_IWSubstring* token : tokens) {
    if (token->empty()) {
      cerr << "CoordinateBox::BuildFromSmilesToken:empty token(s) " << from_smiles << '\n';
      return 0;
    }
  }

  if (! tokens[0]->numeric_value(_cell_size) || _cell_size <= 0.0) {
    cerr << "CoordinateBox::BuildFromSmilesToken:invalid resolution " << from_smiles << "'\n";
    return 0;
  }

  double dims[2];
  for (int i = 0; i < 2; ++i) {
    if (! tokens[i + 1]->numeric_value(dims[i]) || dims[i] <= 0.0) {
      cerr << "CoordinateBox::BuildFromSmilesToken:invalid cell dimensions " << from_smiles << '\n';
      return 0;
    }
  }

  // Make sure that the dimensions are roughly divisible by _cell_size;
  for (int i = 0; i < 2; ++i) {
    double hopefully_integer = dims[i] / _cell_size;
    if (hopefully_integer - static_cast<int>(hopefully_integer + 0.00001) > 0.1) {
      cerr << "CoordinateBox::BuildFromSmilesToken:invalid cell size " << _cell_size << " for dimension " << dims[i] << '\n';
      return 0;
    }
    _cell[i] = static_cast<int>(hopefully_integer);
  }

  return 1;
}

// Each cell_number maps to a unique position in the box. The first task here is
// to identify the 3 positions defined by `cell_number`.
// Note that Z is unconstrained.
template <typename T> 
std::tuple<T, T, T>
CoordinateBox::CoordinatesAsTuple(int cell_number) const {
  assert(cell_number >= 0);

  const int zcell = cell_number / (_cell[0] * _cell[1]);
  cell_number = cell_number % (_cell[0] * _cell[1]);
  const int ycell = cell_number / _cell[0];
  const int xcell = cell_number % _cell[0];

  const T x = _origin.x() + _cell_size * xcell;
  const T y = _origin.y() + _cell_size * ycell;
  const T z = _origin.z() + _cell_size * zcell;
  return std::tuple<T, T, T>{x, y, z};
}
template std::tuple<float, float, float> CoordinateBox::CoordinatesAsTuple(int cell_number) const;

template <typename T>
Space_Vector<T>
CoordinateBox::CoordinatesAsVector(int cell_number) const {
  auto [x, y, z] = CoordinatesAsTuple<T>(cell_number);
  return Space_Vector<T>(x, y, z);
}

template Space_Vector<float> CoordinateBox::CoordinatesAsVector(int cell_number) const;
template Space_Vector<double> CoordinateBox::CoordinatesAsVector(int cell_number) const;

// Note that we do not check for out of range X or Y values.
// Can be enabled via an ifdef
template < typename T>
int
CoordinateBox::CellNumber(T x, T y, T z) const {
  int xcount = static_cast<int>(x / _cell_size + 0.49999);
  int ycount = static_cast<int>(y / _cell_size + 0.49999);
  int zcount = static_cast<int>(z / _cell_size + 0.49999);

#ifdef CELL_NUMBER_CHECK_IN_BOUNDS
  if (x < 0.0 || xcount >= _cell[0])
    return -1;
  if (y < 0.0 || ycount >= _cell[1])
    return -1;
  if (z < 0.0)
    return -1;
#endif

  if (xcount >= _cell[0] || ycount >= _cell[1]) {
    cerr << "CoordinateBox::CellNumber:value out of range " << x << ',' << y << ',' << z << '\n';
    return -1;
  }

  return (_cell[0] * _cell[1]) * zcount +
          _cell[0] * ycount +
          xcount;
}

template int CoordinateBox::CellNumber(float x, float y, float z) const;
template int CoordinateBox::CellNumber(double x, double y, double z) const;
}  // namespace coordinate_box
