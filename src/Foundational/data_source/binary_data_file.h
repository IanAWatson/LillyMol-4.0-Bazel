#ifndef FOUNDATIONAL_DATA_SOURCE_BINARY_DATA_FILE_H
#define FOUNDATIONAL_DATA_SOURCE_BINARY_DATA_FILE_H

#include <optional>

#include "Foundational/iwstring/iwstring.h"

namespace binary_data_file {

class BinaryDataFileReader {
  private:
    int _fd;

    bool _good;

    bool _eof;

    unsigned int _next;

    resizable_array<char> _read_buffer;

    int _items_read;

  // private functions

  void DefaultValues();

  bool OpenFile(const char * fname);

  unsigned int GetNextSize() const;

  bool FillReadBuffer();

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

    int SetBufferSize(int buf_size);

    int items_read() const { return _items_read;}

    std::optional<const_IWSubstring> Next(); 
};

class BinaryDataFileWriter {
  private:
    IWString_and_File_Descriptor _output;

  public:
    BinaryDataFileWriter(int fd);

    int Open(const char * fname);
    int Open(IWString& fname);
    int Open(const const_IWSubstring& fname);

    int Close();

    ssize_t Write(const void * data, int nbytes);
    ssize_t Write(const const_IWSubstring& data);
    ssize_t Write(const IWString& data);

    int good() const {
      return _output.good();
    }
};

void SetDefaultReadBufferSize(int value);
void SetDefaultWriteBufferSize(int value);

}  //  namespace binary_data_file

#endif // FOUNDATIONAL_DATA_SOURCE_BINARY_DATA_FILE_H
