#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "binary_data_file.h"

namespace binary_data_file {

using std::cerr;

unsigned int default_read_buffer_size = 4096; 

// We writ 32 bit item sizes.
constexpr uint32_t sizeof_count = sizeof(uint32_t);

void SetDefaultReadBufferSize(int value) {
  if (value < static_cast<int>(sizeof_count)) {
    value = sizeof_count;
  }
  default_read_buffer_size = value;
}

unsigned int default_write_buffer_size = 8192;
void SetDefaultWriteBufferSize(int value) {
  default_write_buffer_size = value;
}

void
BinaryDataFileReader::DefaultValues() {
  _fd = -1;
  _good = true;
  _eof = false;
  _next = 0;

  _read_buffer.resize(default_read_buffer_size);

  _items_read = 0;
}

BinaryDataFileReader::BinaryDataFileReader() {
  DefaultValues();
}

BinaryDataFileReader::BinaryDataFileReader(const char * fname) {
  DefaultValues();
  OpenFile(fname);
}

BinaryDataFileReader::BinaryDataFileReader(IWString& fname) {
  DefaultValues();
  OpenFile(fname.null_terminated_chars());
}

BinaryDataFileReader::BinaryDataFileReader(const const_IWSubstring& fname) {
  DefaultValues();
  IWString tmp(fname);
  OpenFile(tmp.null_terminated_chars());
}

BinaryDataFileReader::~BinaryDataFileReader() {
}

int
BinaryDataFileReader::Open(const char * fname) {
  return OpenFile(fname);
}

int
BinaryDataFileReader::Open(IWString& fname) {
  return OpenFile(fname.null_terminated_chars());
}

int
BinaryDataFileReader::Open(const const_IWSubstring& fname) {
  IWString tmp(fname);
  return OpenFile(tmp.null_terminated_chars());
}

bool
BinaryDataFileReader::OpenFile(const char * fname) {
  if (_fd >= 0) {
    cerr << "BinaryDataFileReader::OpenFile:already open " << _fd << ", no action\n";
    return 0;
  }

  _fd = IW_FD_OPEN(fname, O_RDONLY);
  if (_fd < 0) {
    _good = false;
    return false;
  }

  _good = true;
  _eof = false;
  return true;
}

// Parse the bytes at _read_buffer + _next as an unsigned int.
unsigned int
BinaryDataFileReader::GetNextSize() const {
  const char * cptr = _read_buffer.rawdata() + _next;
  const unsigned int * iptr = reinterpret_cast<const unsigned int*>(cptr);
#ifdef DEBUG_BINARY_DATA_FILE_READER_NEXT
  cerr << "GetNextSize _next " << _next << " size " << (*iptr) << '\n';
#endif
  return *iptr;
}

std::optional<const_IWSubstring>
BinaryDataFileReader::Next() {
  if (_fd < 0 || _eof || ! _good) {
    return std::nullopt;
  }

  // First task is to read the size of the next item.
  // At a minimum, we need 4 bytes for size of next. We
  // must not assume next has non zero size...
  if(_read_buffer.size() - _next < sizeof_count) {
    if (! FillReadBuffer()) {
      return std::nullopt;
    }
    if (_read_buffer.size() - _next < sizeof_count) {
      _good = 0;
      return std::nullopt;
    }
  }

  const int bytes_in_next = GetNextSize();
  _next += sizeof_count;
#ifdef DEBUG_BINARY_DATA_FILE_READER_NEXT
  cerr << "After GetNextSize update _next " << _next << " next item size " << bytes_in_next << '\n';
#endif
  if (bytes_in_next == 0) {
    return const_IWSubstring();
  }

  // If the next object is wholly contained in the read buffer,
  // use it.
  if (_next + bytes_in_next <= _read_buffer.size()) {
    const_IWSubstring result(_read_buffer.rawdata() + _next, bytes_in_next);
    _next += bytes_in_next;
    return result;
  }

  // Will need to read more data, remove all the stuff we have previously consumed.
  _read_buffer.remove_from_to(0, _next);
  _next = 0;

  if (_read_buffer.elements_allocated() < bytes_in_next) {
    _read_buffer.resize(bytes_in_next + sizeof_count);
  }

  if (! FillReadBuffer()) {
    _good = 0;
    return std::nullopt;
  }

  _next = bytes_in_next;

  const_IWSubstring result(_read_buffer.rawdata(), bytes_in_next);
  return result;
}

bool
BinaryDataFileReader::FillReadBuffer() {
  const int existing_size = _read_buffer.number_elements();

  unsigned int unused_capacity = _read_buffer.elements_allocated() - existing_size;
  if (unused_capacity < default_read_buffer_size) {
    _read_buffer.resize(_read_buffer.elements_allocated() + default_read_buffer_size);
    unused_capacity = _read_buffer.elements_allocated() - existing_size;
  }
  int bytes_read = IW_FD_READ(_fd, _read_buffer.rawdata() + existing_size, unused_capacity);
  if (bytes_read < 0) {
    _good = false;
    return false;
  }
  if (bytes_read == 0) {
    _eof = true;
    return false;
  }

  _read_buffer.set_number_elements(existing_size + bytes_read);
  return true;
}

BinaryDataFileWriter::BinaryDataFileWriter(int fd) : _output(fd) {
}

int
BinaryDataFileWriter::Open(const char * fname) {
  if (_output.is_open()) {
    std::cerr << "BinaryDataFileWriter::Open:already open\n";
    return 0;
  }
  return _output.open(fname);
}

int
BinaryDataFileWriter::Open(IWString& fname) {
  if (_output.is_open()) {
    std::cerr << "BinaryDataFileWriter::Open:already open\n";
    return 0;
  }
  return _output.open(fname.null_terminated_chars());
}

int
BinaryDataFileWriter::Open(const const_IWSubstring& fname) {
  if (_output.is_open()) {
    std::cerr << "BinaryDataFileWriter::Open:already open\n";
    return 0;
  }
  IWString tmp(fname);
  return _output.open(tmp.null_terminated_chars());
}

ssize_t
BinaryDataFileWriter::Write(const void * data, int nbytes) {
//cerr << "BinaryDataFileWriter::Write: writing " << nbytes << " bytes\n";
  const uint32_t mysize = nbytes;
  const char * cptr = reinterpret_cast<const char *>(&mysize);
  if (_output.write(cptr, sizeof_count) != sizeof_count) {
    std::cerr << "BinaryDataFileWriter::Write:size not written\n";
    return 0;
  }
  const auto bytes_written = _output.write(reinterpret_cast<const char *>(data), nbytes);
  if (bytes_written != nbytes) {
    std::cerr << "BinaryDataFileWriter::Write:wrote " << bytes_written << " of " << nbytes << " bytes\n";
    return 0;
  }
  _output.write_if_buffer_holds_more_than(default_write_buffer_size);
  return sizeof(unsigned int) + bytes_written;
}

ssize_t
BinaryDataFileWriter::Write(const IWString& data) {
  return Write(data.rawdata(), data.size());
}
ssize_t
BinaryDataFileWriter::Write(const const_IWSubstring& data) {
  return Write(data.data(), data.length());
}

int
BinaryDataFileWriter::Close() {
  return _output.close();
}

}  //  namespace binary_data_file
