#include <functional>

#include "misc.h"

namespace iwmisc {
RestoreState::RestoreState(getter g, setter s, int state) : _setter(s) {
  _initial_state = std::invoke(g);
  std::invoke(s, state);
}

RestoreState::~RestoreState() {
  std::invoke(_setter, _initial_state);
}

}  // namespace iwmisc
