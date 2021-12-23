// Read and write TFDataRecord

#include <optional>

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/iwstring.h"

namespace iw_tf_data_record {
class TFDataWriter {
  private:
  public:
};

class TFDataReader {
  private:
    // The file descriptor from which data is retrieved.
    int _fd;
    // Status flags.

    bool _good;
    bool _eof;

    int _items_read;

    // An index into _read_buffer where the next item starts.
    uint64_t _next;

    // Data is read from _fd into _read_buffer.
    // Note that the resizable_array class is not 64 bit capable.
    resizable_array<char> _read_buffer;

    // private functions.

    void DefaultValues();

    bool OpenFile(const char * fname);

    std::optional<uint64_t> GetLength();

    bool FillReadBuffer(uint64_t bytes_needed = 8192);

  public:
    TFDataReader();
    TFDataReader(const char * fname);
    TFDataReader(IWString & fname);
    TFDataReader(const const_IWSubstring & fname);
    ~TFDataReader();

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
    // the next item.
    std::optional<const_IWSubstring> Next(); 

    // Read serialized proto of type P and return decoded form.
    template <typename P>
    std::optional<P>
    ReadProto();
};
}  // namespace iw_tf_data_record
