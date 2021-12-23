#include <algorithm>
#include <cstdint>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "crc32c/crc32c.h"

#include "tfdatarecord.h"

namespace iw_tf_data_record {

using std::cerr;

unsigned int default_read_buffer_size = 4096; 

constexpr uint64_t sizeof_length = sizeof(uint64_t);
constexpr uint64_t sizeof_crc = sizeof(uint32_t);

void
TFDataReader::DefaultValues() {
  _fd = -1;
  _good = true;
  _eof = false;
  _next = 0;

  _read_buffer.resize(default_read_buffer_size);

//_compression_type = kUncompressed;

  _items_read = 0;
}

TFDataReader::TFDataReader() {
  DefaultValues();
}

TFDataReader::TFDataReader(const char * fname) {
  DefaultValues();
  OpenFile(fname);
}

TFDataReader::TFDataReader(IWString& fname) {
  DefaultValues();
  OpenFile(fname.null_terminated_chars());
}

TFDataReader::TFDataReader(const const_IWSubstring& fname) {
  DefaultValues();
  IWString tmp(fname);
  OpenFile(tmp.null_terminated_chars());
}

TFDataReader::~TFDataReader() {
}

int
TFDataReader::Open(const char * fname) {
  return OpenFile(fname);
}

int
TFDataReader::Open(IWString& fname) {
  return OpenFile(fname);
}

int
TFDataReader::Open(const const_IWSubstring& fname) {
  IWString tmp(fname);
  return OpenFile(tmp);
}

bool
TFDataReader::OpenFile(const char * fname) {
  if (_fd >= 0) {
    cerr << "TFDataReader::OpenFile:already open " << _fd << ", no action\n";
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

std::optional<uint64_t>
TFDataReader::GetLength() {
  const uint64_t* lptr = reinterpret_cast<const uint64_t*>(_read_buffer.rawdata() + _next);
  const uint32_t* crc = reinterpret_cast<const uint32_t*>(_read_buffer.rawdata() + _next + sizeof_length);

  uint32_t result = crc32c::Crc32c(reinterpret_cast<const uint8_t*>(lptr), sizeof_length);
  if (result != *crc) {
    cerr << "TFDataReader::GetLength:crc fails, length " << *lptr << '\n';
    _good = 0;
    return std::nullopt;
  }
  _next += sizeof_length + sizeof_crc;
  return *lptr;
}

std::optional<const_IWSubstring>
TFDataReader::Next() {
  if (_eof || ! _good) {
    cerr << "TFDataReader::Next:eof or not good\n";
    return std::nullopt;
  }

  // First task is to read the size of the next item.
  // At a minimum, we need 8+4 bytes for size of next and the crc
  if(_read_buffer.size() - _next < (sizeof_length + sizeof_crc)) {
    if (! FillReadBuffer()) {
      return std::nullopt;
    }
    if (_read_buffer.size() - _next < (sizeof_length + sizeof_crc)) {
      _good = 0;
      return std::nullopt;
    }
  }

  std::optional<uint64_t> length = GetLength();
  if (! length) {
    return std::nullopt;
  }
  if (_read_buffer.size() - _next + *length + sizeof_crc > _read_buffer.size()) {
    if (! FillReadBuffer(*length)) {
      return std::nullopt;
    }
  }

  const_IWSubstring result(_read_buffer.rawdata() + _next, *length);

  const uint32_t crc_data = crc32c::Crc32c(reinterpret_cast<const uint8_t*>(result.data()), *length);
  const uint32_t * crc = reinterpret_cast<const uint32_t*>(_read_buffer.rawdata() + _next + *length);
  if (crc_data != *crc) {
    cerr << "TFDataReader::Next:Invalid data crc " << *length << " bytes\n";
    _good = 0;
    return std::nullopt;
  }

  ++_items_read;
  return result;
}

#ifdef NOT_NEEDED_JJJ
bool
TFDataReader::FillReadBuffer() {
  const int existing_size = _read_buffer.number_elements();

//cerr << "FillReadBuffer: starting pos " << IW_FD_LSEEK(_fd, 0, SEEK_CUR) << " bytes\n";
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
#endif

// The class is needs to be able to read `bytes_needed` into an item.
// Upon exit, the next bytes_needed of data will be in the buffer.
bool
TFDataReader::FillReadBuffer(uint64_t bytes_needed) {
  // If the next bytes_needed bytes are already in _read_buffer we are done.
  if (_next + bytes_needed <= _read_buffer.number_elements()) {
    return true;
  }
  
  // We need to shift the data and maybe resize.
  _read_buffer.remove_from_to(0, _next);
  _next = 0;
  if (_read_buffer.elements_allocated() < bytes_needed) {
    _read_buffer.resize(bytes_needed);
  }

  int to_read = std::max(default_read_buffer_size, static_cast<uint32_t>(_read_buffer.elements_allocated()));

  int bytes_read = IW_FD_READ(_fd, _read_buffer.rawdata(), to_read);
  if (bytes_read < 0) {
    _good = false;
    return false;
  }
  if (bytes_read == 0) {
    _eof = true;
    return false;
  }
  if (bytes_read != to_read) {
    cerr << "TFDataReader::cannot read " << to_read << " bytes\n";
    _good = false;
    return false;
  }

  _read_buffer.set_number_elements(bytes_read);
  return true;
}


}  // namespace iw_tf_data_record
