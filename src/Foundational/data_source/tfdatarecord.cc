#include <algorithm>
#include <cstdint>

#include <endian.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "crc32c/crc32c.h"

#include "Foundational/iwmisc/misc.h"

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
#ifdef DEBUG_IW_TF_DATA
  cerr << "Openfile " << fname << " returned " << _fd << '\n';
#endif
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

#ifdef DEBUG_IW_TF_DATA
  cerr << "Length " << *lptr  << " crc " << *crc << '\n';
#endif
  uint64_t len = *lptr;
  uint32_t result = crc32c::Crc32c(reinterpret_cast<const char*>(&len), sizeof_length);
  result = ((result >> 15) | (result << 17)) + 0xa282ead8ul;
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
      _eof = 1;
      return std::nullopt;
    }
  }

  std::optional<uint64_t> length = GetLength();
  if (! length) {
    return std::nullopt;
  }
  if (_next + *length + sizeof_crc > static_cast<uint64_t>(_read_buffer.number_elements())) {
    if (! FillReadBuffer(*length + sizeof_crc)) {
      return std::nullopt;
    }
  }

  const_IWSubstring result(_read_buffer.rawdata() + _next, *length);

  const uint32_t crc_data = crc32c::Crc32c(result.data(), *length);
  const uint32_t masked_crc = ((crc_data >> 15) | (crc_data << 17)) + 0xa282ead8ul;
  const uint32_t * crc = reinterpret_cast<const uint32_t*>(_read_buffer.rawdata() + _next + *length);
  if (*crc != masked_crc) {
    cerr << "TFDataReader::Next:Invalid data crc " << *length << " bytes\n";
    _good = 0;
    return std::nullopt;
  }

  _next += *length + sizeof_crc;

  ++_items_read;
  return result;
}

// The class is needs to be able to read `bytes_needed` into an item.
// Upon exit, the next bytes_needed of data will be in the buffer.
bool
TFDataReader::FillReadBuffer(uint64_t bytes_needed) {
  // If the next bytes_needed bytes are already in _read_buffer we are done.
  if (_next + bytes_needed <= static_cast<uint64_t>(_read_buffer.number_elements())) {
    return true;
  }
  
  // We need to shift the data and maybe resize.
  if (_next > 0) {
    _read_buffer.remove_from_to(0, _next);
    _next = 0;
  }
  const int bytes_already_present = _read_buffer.number_elements();

  int to_read = default_read_buffer_size;
  if (static_cast<uint64_t>(_read_buffer.elements_allocated()) < bytes_needed) {
    _read_buffer.resize(bytes_needed);
    to_read = bytes_needed - _read_buffer.number_elements();
  }

  int bytes_read = IW_FD_READ(_fd, _read_buffer.rawdata() + bytes_already_present, to_read);
#ifdef DEBUG_IW_TF_DATA
  cerr << "Reading " << to_read << " got " << bytes_read << " bytes from file " << _fd << '\n';
#endif
  if (bytes_read < 0) {
    _good = false;
    return false;
  }
  if (bytes_read == 0) {
    _eof = true;
    return false;
  }
  if (bytes_read != to_read) {
    _eof = true;
  }

  _read_buffer.set_number_elements(bytes_already_present + bytes_read);
  return true;
}


}  // namespace iw_tf_data_record
