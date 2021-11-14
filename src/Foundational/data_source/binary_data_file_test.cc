// Tester for binary data file

#include <sys/stat.h>
#include <fcntl.h>

#include <cstdio>
#include <random>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "binary_data_file.h"

namespace {
// If we ever change from 32 bit item sizes, this will
// need to change. It could possibly be shared with the
// class source file, but it is an internal detail.
constexpr uint32_t sizeof_count = sizeof(uint32_t);

using std::cerr;

using binary_data_file::BinaryDataFileReader;
using binary_data_file::BinaryDataFileWriter;

class TestBinaryDataFileReader : public testing::Test {
  protected:
    IWString _fname;
    int _fd;

    void SetUp();
    void TearDown();
};

void
TestBinaryDataFileReader::SetUp() {
  _fname << "/tmp/bdftestXXXXXX.dat";
  _fd = ::mkostemps(const_cast<char*>(_fname.null_terminated_chars()), 4, O_WRONLY | O_CREAT | O_TRUNC);
  if (_fd < 0) {
    cerr << "TestBinaryDataFileReader::SetUp:mkostemps failed, rc " << _fd << '\n';
  }
}

void
TestBinaryDataFileReader::TearDown() {
  ::remove(_fname.null_terminated_chars());
}

resizable_array<int>
ArrayOfSize(int s) {
  resizable_array<int> result(s);
  for (int i = 0; i < s; ++i) {
    result << i;
  }
  return result;
}

TEST_F(TestBinaryDataFileReader, TestSingleItem) {
  BinaryDataFileWriter writer(_fd);
  constexpr int n = 1;
  const resizable_array<int> data = ArrayOfSize(n);
  int bytes_written = writer.Write(data.rawdata(), n * sizeof(int));
  ASSERT_EQ(bytes_written, sizeof_count + sizeof(int));
  writer.Close();

  binary_data_file::SetDefaultReadBufferSize(2);
  binary_data_file::BinaryDataFileReader bdf(_fname);
  ASSERT_TRUE(bdf.good());
  std::optional<const_IWSubstring> as_read = bdf.Next();
  ASSERT_TRUE(as_read);
  ASSERT_FALSE(as_read->empty());
  ASSERT_EQ(as_read->length(), sizeof(int));
  const int * as_int = reinterpret_cast<const int *>(as_read->data());
  EXPECT_EQ(*as_int, 0);
}

TEST_F(TestBinaryDataFileReader, TestSingleItemSizeTwo) {
  BinaryDataFileWriter writer(_fd);
  constexpr int n = 2;
  const resizable_array<int> data = ArrayOfSize(n);
  int bytes_written = writer.Write(data.rawdata(), n * sizeof(int));
  ASSERT_EQ(bytes_written, sizeof_count + n * sizeof(int));
  writer.Close();

  binary_data_file::SetDefaultReadBufferSize(2);
  binary_data_file::BinaryDataFileReader bdf(_fname);
  ASSERT_TRUE(bdf.good());
  std::optional<const_IWSubstring> as_read = bdf.Next();
  ASSERT_TRUE(as_read);
  ASSERT_FALSE(as_read->empty());
  ASSERT_EQ(as_read->length(), 2 * sizeof(int));
  const int * as_int = reinterpret_cast<const int *>(as_read->data());
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(as_int[i], i);
  }
}

TEST_F(TestBinaryDataFileReader, TestTwoItems1) {
  BinaryDataFileWriter writer(_fd);
  constexpr int nitems = 2;
  const resizable_array<int> data = ArrayOfSize(1);
  ASSERT_EQ(writer.Write(data.rawdata(), data.size() * sizeof(int)), sizeof_count + 1 * sizeof(int));
  ASSERT_EQ(writer.Write(data.rawdata(), data.size() * sizeof(int)), sizeof_count + 1 * sizeof(int));
  writer.Close();

  BinaryDataFileReader reader(_fname);

  int items_retrieved = 0;
  do {
    auto item = reader.Next();
    if (! item) {
      break;
    }
    ++items_retrieved;
  } while (1);

  EXPECT_EQ(items_retrieved, nitems);
}

TEST_F(TestBinaryDataFileReader, TestLots) {
  constexpr int n = 3000;
  BinaryDataFileWriter writer(_fd);

  resizable_array<int> sizes(n);
  for (int i = 0; i < n; ++i) {
    sizes << i;
  }
  std::random_shuffle(sizes.begin(), sizes.end());
  // Make sure the last item in there is zero size.
  sizes << 0;
  // And add a couple of consecutive zero size things
  sizes[50] = 0;
  sizes[51] = 0;
  sizes[52] = 0;

  for (int s : sizes) {
    const resizable_array<int> data = ArrayOfSize(s);
    ssize_t expected = sizeof_count + s * sizeof(int);
    ASSERT_EQ(writer.Write(data.rawdata(), s * sizeof(int)), expected);
  }
  writer.Close();

  binary_data_file::SetDefaultReadBufferSize(6);
  BinaryDataFileReader reader(_fname);
  ASSERT_TRUE(reader.good());

  int nitems = sizes.number_elements();
  for (int i = 0; i < nitems; ++i) {
    auto item = reader.Next();
    ASSERT_TRUE(item);
    ASSERT_EQ(item->length(), sizes[i] * sizeof(int));
    const int * as_int = reinterpret_cast<const int *>(item->data());
    for (int j = 0; j < sizes[i]; ++j) {
      EXPECT_EQ(as_int[j], j);
    }
  }
  EXPECT_FALSE(reader.eof());
  EXPECT_FALSE(reader.Next());
  EXPECT_TRUE(reader.eof());
  EXPECT_FALSE(reader.Next());
  EXPECT_TRUE(reader.eof());
}

}  // namespace
