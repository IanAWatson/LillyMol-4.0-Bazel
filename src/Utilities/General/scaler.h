#ifndef UTILITIES_GENERAL_SCALER_H_
#define UTILITIES_GENERAL_SCALER_H_

// Rescaling floating point values to and from the [0,1] range.

#include "Utilities/General/feature_scaling.pb.h"

namespace feature_scaler {

template <typename T>
class FeatureScaler {
  private:
    T _min;
    T _max;
    T _range;

    // What to do if a value is outside the range.
    bool _truncate_out_of_range;

  // private functions

    void _default_values();

  public:
    FeatureScaler();

    // Initialized with a min and max.
    FeatureScaler(T zmin, T zmax);

    int Initialise(const FeatureScaling::FeatureScaling& proto);

    int Active() const {
      return _range > T();
    }

    void set_truncate_out_of_range(bool s) {
      _truncate_out_of_range = s;
    }

    // Factory construction from a proto.
    static FeatureScaler Build(const FeatureScaling::FeatureScaling& proto);

    bool InRange(T value) const {
      return value >= _min && value <= _max;
    }

    // Once built, the scaling operations supported.
    T ScaleTo01(T value) const;
    T ScaleBackToOrignalRange(T value) const;
};

#ifdef FEATURE_SCALER_IMPLEMENTATION

template <typename T>
void
FeatureScaler<T>::_default_values() {
  _range = T();
  _truncate_out_of_range = false;
}

template <typename T>
FeatureScaler<T>::FeatureScaler() {
  _default_values();

  _min = T();
  _max = T();
  _range = T();
}

template <typename T>
FeatureScaler<T>::FeatureScaler(T zmin, T zmax) {
  _default_values();

  _min = zmin;
  _max = zmax;
  _range = _max - _min;
}

template <typename T>
FeatureScaler<T>
FeatureScaler<T>::Build(const FeatureScaling::FeatureScaling& proto) {
  FeatureScaler<T> result;
  result._default_values();

  result._min = proto.min();
  result._max = proto.max();
  result._range = result._max - result._min;

  return result;
}

template <typename T>
int
FeatureScaler<T>::Initialise(const FeatureScaling::FeatureScaling& proto) {
  _default_values();
  _min = proto.min();
  _max = proto.max();
  _range = _max - _min;
  return 1;
}

template <typename T>
T
FeatureScaler<T>::ScaleTo01(T value) const {
  if (value >= _min && value <= _max) {
    return (value - _min) / _range;
  }

#ifdef STATS_ON_OUT_OF_RANGE
  // Out of range, update stats and decide what to do.
  if (value < scaling.min()) {
    out_of_range_low++;
  } else {
    out_of_range_high++;
  }
#endif

  if (_truncate_out_of_range) {
    if (value < _min) {
      return T();
    } else {
      return 1;
    }
  }

  // Write extrapolated value.
  if (value < _min) {
    return (value - _min) / _range;
  } else {
    return 1.0 + (value - _max) / _range;
  }
}

template <typename T>
T
FeatureScaler<T>::ScaleBackToOrignalRange(T value) const {
  assert(value >= T() && value <= 1.0);

  if (value >= T() && value <= 1.0) {
    return _min + value * _range;
  }

#ifdef STATS_ON_OUT_OF_RANGE
  // Out of range, update stats and decide what to do.
  if (value < scaling.min()) {
    out_of_range_low++;
  } else {
    out_of_range_high++;
  }
#endif

  if (_truncate_out_of_range) {
    if (value < T()) {
      return _min;
    } else {
      return _max;
    }
  }

  // Write extrapolated value.
  if (value < T()) {
    return _min + value * _range;
  } else {
    return _max + (value - 1.0) * _range;
  }
}

#endif   // FEATURE_SCALER_IMPLEMENTATION

}  // namespace feature_scaler

#endif  // UTILITIES_GENERAL_SCALER_H_
