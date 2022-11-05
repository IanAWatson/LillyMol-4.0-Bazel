// AFile implementation functions.

#include <unistd.h>

#include "Foundational/iwstring/iwstring.h"

#include "proto_support.h"

namespace iwmisc {

AFile::AFile(IWString& fname, int mode) {
  int flags = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  _fd = IW_FD_OPEN(fname.null_terminated_chars(), mode, flags);
}

AFile::~AFile() {
  IW_FD_CLOSE(_fd);
}

size_t
AFile::ReadAll(IWString& buffer) {
  const off_t current = lseek(_fd, 0, SEEK_CUR);
  const off_t remaining = lseek(_fd, 0, SEEK_END);
  lseek(_fd, current, SEEK_SET);

  buffer.extend(remaining);

  if (auto x =  read(_fd, (void*) buffer.data(), remaining); x != remaining) {
    cerr << "AFile::ReadAll:cannot read " << remaining << " bytes, got " << x << " err " << errno << '\n';
    return 0;
  }


  return remaining;
}


}  // namespace iwmisc
