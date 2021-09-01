#ifndef FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
#define FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
// Functions to support operating with protos.

#include <ostream>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/iwstring/iwstring.h"

namespace iwmisc {
// Write a proto to `fname` using Text_Format.
template <typename Proto>
int
WriteProtoAsText(const Proto& proto, IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    std::cerr << "WriteProtoAsText:cannot open " << fname << '\n';
    return 0;
  }

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(proto, zero_copy_output.get())) {
    std::cerr << "WriteProtoAsText:cannot write " << fname << "\n";
    return 0;
  }

  return 1;
}

}  // namespace iwmisc

#endif // FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
