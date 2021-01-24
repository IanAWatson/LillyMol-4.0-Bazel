#include "Foundational/iwmisc/iwre2.h"

namespace iwre2 {
bool RE2FullMatch(const IWString& s, RE2& rx) {
  re2::StringPiece tmp(s.data(), s.length());
  return RE2::FullMatch(tmp, rx);
}
bool RE2FullMatch(const const_IWSubstring& s, RE2& rx) {
  re2::StringPiece tmp(s.data(), s.length());
  return RE2::FullMatch(tmp, rx);
}
bool RE2PartialMatch(const IWString& s, RE2& rx) {
  re2::StringPiece tmp(s.data(), s.length());
  return RE2::PartialMatch(tmp, rx);
}
bool RE2PartialMatch(const const_IWSubstring& s, RE2& rx) {
  re2::StringPiece tmp(s.data(), s.length());
  return RE2::PartialMatch(tmp, rx);
}

bool RE2Reset(std::unique_ptr<RE2>& rx, const IWString& s) {
  re2::StringPiece tmp(s.data(), s.length());
  rx.reset(new RE2(tmp));
  return rx->ok();
}
bool RE2Reset(std::unique_ptr<RE2>& rx, const const_IWSubstring& s) {
  re2::StringPiece tmp(s.data(), s.length());
  rx.reset(new RE2(tmp));
  return rx->ok();
}
}  // namespace iwre2
