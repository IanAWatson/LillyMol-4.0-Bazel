#ifndef FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
#define FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
// Functions to support operating with protos.

#include <fcntl.h>
#include <optional>
#include <ostream>

#include "absl/strings/string_view.h"

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/stubs/stringpiece.h"
#include "google/protobuf/util/json_util.h"

#include "Foundational/iwstring/iwstring.h"

namespace iwmisc {

using std::cerr;

// A lightweight class to open a file descrptor and make sure
// it gets closed.
class AFile {
  private:
    int _fd;

  public:
    AFile(IWString& fname, int mode);  // O_RDONLY or O_WRONLY
    ~AFile();

    int good() const {
      return _fd >= 0;
    }

    int fd() const {
      return _fd;
    }

    size_t ReadAll(IWString& buffer);
};

// Write a proto to `fname` using Text_Format.
template <typename Proto>
int
WriteProtoAsText(const Proto& proto, IWString& fname) {
  AFile output(fname, O_WRONLY | O_TRUNC | O_CREAT);
  if (! output.good()) {
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

template <typename Proto>
int
WriteBinaryProto(const Proto& proto, IWString& fname) {
  AFile output(fname, O_WRONLY | O_TRUNC | O_CREAT);
  if (! output.good()) {
    cerr << "WriteBinaryProto:cannot open '" << fname << "'\n";
    return 0;
  }

  return proto.SerializeToFileDescriptor(output.fd());
}

template <typename Proto>
std::optional<Proto>
ReadBinaryProto(IWString& fname) {
  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadBinaryProto:cannot open '" << fname << "'\n";
    return std::nullopt;
  }

  Proto result;
  if (! result.ParseFromFileDescriptor(input.fd())) {
    cerr << "ReadBinaryProto:cannot parse '" << fname << "'\n";
    return std::nullopt;
  }

  return result;
}

// Read existing proto.
template <typename Proto>
bool
ReadBinaryProto(IWString& fname, Proto& proto) {
  proto.Clear();

  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadBinaryProto:cannot open '" << fname << "'\n";
    return false;
  }

  if (! proto.ParseFromFileDescriptor(input.fd())) {
    cerr << "ReadBinaryProto:cannot parse '" << fname << "'\n";
    return false;
  }

  return true;
}

template <typename Proto>
std::optional<Proto>
ReadTextProto(IWString& fname) {
  AFile input(fname, O_RDONLY);
  if (! input.good()) {
    cerr << "ReadTextProto:cannot open '" << fname << "'\n";
    return std::nullopt;
  }

  using google::protobuf::io::ZeroCopyInputStream;
  using google::protobuf::io::FileInputStream;
  std::unique_ptr<FileInputStream> zero_copy_input(new FileInputStream(input.fd()));

  Proto result;

  if (fname.ends_with(".json")) {
    IWString data;
    if (!input.ReadAll(data)) {
      cerr << "ReadTextProto:cannot read\n";
      return std::nullopt;
    }
    //google::protobuf::stringpiece_internal::StringPiece string_piece(data.data(), data.size());
    absl::string_view string_piece(data.data(), data.size());
    google::protobuf::util::JsonParseOptions options;
    auto status = google::protobuf::util::JsonStringToMessage(string_piece, &result, options);
    if (! status.ok()) {
      cerr << "ReadTextProto:cannot read json '" << fname << " " << status << '\n';
      return std::nullopt;
    }
  } else {
    if (! google::protobuf::TextFormat::Parse(zero_copy_input.get(), &result)) {
      cerr << "ReadTextProto:cannot read '" << fname << "'\n";
      return std::nullopt;
    }
  }

  return result;
}

template <typename Proto>
int
WriteTextProto(Proto& proto, IWString& fname) {
  AFile output(fname, O_WRONLY | O_TRUNC | O_CREAT);

  if (! output.good()) {
    cerr << "WriteTextProto:cannot open '" << fname << "'\n";
    return 0;
  }

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(proto, zero_copy_output.get())) {
    cerr << "WriteTextProto:cannot write '" << fname << "'\n";
    return 0;
  }

  return 1;
}

}  // namespace iwmisc

#endif // FOUNDATIONAL_IWMISC_PROTO_SUPPORT_H_
