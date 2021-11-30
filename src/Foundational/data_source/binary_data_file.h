#ifndef FOUNDATIONAL_DATA_SOURCE_BINARY_DATA_FILE_H
#define FOUNDATIONAL_DATA_SOURCE_BINARY_DATA_FILE_H

#include <optional>
#include <string>

#include "snappy.h"

#include "Foundational/iwstring/iwstring.h"

namespace binary_data_file {

// Data can be optionally compressed.
enum CompressionType {
  kUncompressed,
  kSnappy,
  // Last entry which is used for determining if a value is valid or not.
  // Insert new compression types before this.
  kInvalid,
};

// The BinaryDataFileWriter and BinaryDataFileReader classes
// write files containing blobs.
// As a blob is stored, the size of the blob is written to the file first.
// This makes it useful for storing things like serialised protos.
// I did look at the library in protobuf that does this, but it seemed
// overly complex for my anticipated use.
// Have not tested this for endianness problems, but hopefully OK.
// Efficiency is the objective, and there should be zero copying of data.
class BinaryDataFileReader {
  private:
    // The file descriptor from which data is retrieved.
    int _fd;

    // Status flags.
    bool _good;
    bool _eof;

    // An index into _read_buffer where the next item starts.
    unsigned int _next;

    // Data is read from _fd into _read_buffer.
    resizable_array<char> _read_buffer;

    // Optional compression.

    CompressionType _compression_type;

    // If compression is active, we need somewhere to hold
    // the uncompressed data, since we always return a const_IWSubstring.
    std::string _uncompressed;

    // For monitoring.
    int _items_read;

  // private functions

  void DefaultValues();

  bool OpenFile(const char * fname);

  unsigned int GetNextSize() const;

  bool FillReadBuffer();

  // The file begins with a magic number (4 bytes). Read that
  // and make sure it is correct.
  // Returns true if the word can be read and is correct.
  bool ReadMagicNumber();
  // Following the header is a 4 byte word for the compression type.
  bool ReadCompressionType();
  // Reads both the magic number and the compression type.
  int ReadHeader();

  // Raw data has been read. If no compression is active, just return
  // `data`. If uncompression is needed, then do that.
  // 
  const_IWSubstring MaybeExpand(const const_IWSubstring& data);

  public:
    BinaryDataFileReader();
    BinaryDataFileReader(const char * fname);
    BinaryDataFileReader(IWString & fname);
    BinaryDataFileReader(const const_IWSubstring & fname);
    ~BinaryDataFileReader();

    int Open(const char * fname);
    int Open(IWString & fname);
    int Open(const const_IWSubstring & fname);

    bool IsOpen() const { return _fd >= 0;}

    bool good() const { return _good;}
    bool eof() const { return _eof;}

    int Close();

    // Will not reduce the size of _read_buffer, but may
    // increase it.
    int SetBufferSize(int buf_size);

    int items_read() const { return _items_read;}

    // The primary method for this class. If possible, return
    // a pointer to the next item.
    std::optional<const_IWSubstring> Next(); 

    // Read serialized proto of type P and return decoded form.
    template <typename P>
    std::optional<P>
    ReadProto();
};

// Writes data files that can subsequently be read by BinaryDataFileReader.
class BinaryDataFileWriter {
  private:
    // Handy abstraction that already has a file descriptor and
    // methods for writing.
    IWString_and_File_Descriptor _output;

    // Optional compression.
    CompressionType _compression_type;

    // And if compressed, a buffer to hold uncompressed data
    std::string _compressed;

  // private functions
    void DefaultValues();

    int WriteMagicNumber();
    int WriteCompressionType();
    // Function that writes both the magic number and the compression type.
    int WriteHeader();

    const_IWSubstring MaybeCompressed(const void * data, int nbytes);

  public:
    BinaryDataFileWriter();
    BinaryDataFileWriter(int fd, CompressionType ctype = kUncompressed);

    int good() const {
      return _output.good();
    }

    // Compression must be specified before the file is opened.
    bool SetCompression(CompressionType compression_type);

    int Open(const char * fname);
    int Open(IWString& fname);
    int Open(const const_IWSubstring& fname);

    int Close();

    ssize_t Write(const void * data, int nbytes);
    ssize_t Write(const const_IWSubstring& data);
    ssize_t Write(const IWString& data);

    // Write a serialized proto.
    template <typename P>
    ssize_t WriteSerializedProto(const P& proto);
};

// Upon construction, the reader and writer classes initialise
// their internal buffers with these values.
void SetDefaultReadBufferSize(int value);
void SetDefaultWriteBufferSize(int value);

// only include if needed.
#ifdef BINARY_FILE_PROTO_IMPLEMENTATION

#include "google/protobuf/message_lite.h"

template <typename P>
ssize_t
BinaryDataFileWriter::WriteSerializedProto(const P& proto) {
  std::string data;
  proto.SerializeToString(&data);
  return Write(data.data(), data.size());
}

template <typename P>
std::optional<P>
BinaryDataFileReader::ReadProto() {
  std::optional<const_IWSubstring> maybe_data = Next();
  if (! maybe_data) {
    return std::nullopt;
  }

  P proto;
  if (! proto.ParseFromArray(maybe_data->data(), maybe_data->length())) {
    std::cerr << "BinaryDataFileReader::ReadProto:cannot parse\n";
    return std::nullopt;
  }

  return proto;
}

#endif // BINARY_FILE_PROTO_IMPLEMENTATION

}  //  namespace binary_data_file

#endif // FOUNDATIONAL_DATA_SOURCE_BINARY_DATA_FILE_H
