# Python tools to read and write files used by
# BinaryDataFileReader and BinaryDataFileWriter

import ctypes
from enum import IntEnum
import snappy

from typing import Optional

MAGIC_NUMBER = 3177520567

class CompressionType(IntEnum):
  ctype_uncompressed = 0
  ctype_snappy = 1
  ctype_invalid = 2

class BinaryDataFileWriter:
  """Create files compatible with Foundational/datasource/binary_data_file.h
  """
  def __init__(self):
    self._stream = None
    self._compression_type = CompressionType.ctype_uncompressed

  def __del__(self):
    self.close()

  def set_compression(self, ctype: CompressionType):
    self._compression_type = ctype

  def open(self, fname: str) -> bool:
    """Open `fname` for writing.

    Args:
      fname: the file to open
    Returns:
      True on success
    """
    self._stream = open(fname, "wb")
    return self.write_number(MAGIC_NUMBER) and self.write_number(int(self._compression_type));

  def close(self):
    self._stream.close()

  def flush(self):
    self._stream.flush();

  def write(self, data):
    if len(data) == 0:
      self.write_number(len(data))
      return True

    if self._compression_type == CompressionType.ctype_uncompressed:
      self.write_number(len(data))
      self._stream.write(str.encode(data))
      return True

    if self._compression_type == CompressionType.ctype_snappy:
      compressed = snappy.compress(data)
      self.write_number(len(compressed))
      self._stream.write(compressed)

    return True
      
  def write_number(self, value: int) -> bool:
    """Write the uint32_t binary value of `i`.
    Args:
      value: to be written
    """
#   print(f"Writing binary {value}")
    uint = ctypes.c_uint(value)
#   print(f"valye {value} written as {uint}")
    self._stream.write(uint)
#   print("Writing binary value returning True")
    return True

class BinaryDataFileReader:
  """Class to read files created by BinaryDataFileWriter.
  """
  def __init__(self, fname: str):
    """
    Args:
      fname: file to read.
    """
    self._compression_type = CompressionType.ctype_uncompressed
    self._stream = open(fname, "rb")

    value = self._stream.read(4)
    value = int.from_bytes(value, "little", signed=False)
    if value != MAGIC_NUMBER:
      raise ValueError(f"Invalid magic number got {value} not MAGIC_NUMBER")

    value = self._stream.read(4)
    value = int.from_bytes(value, "little", signed=False)
    if not value in CompressionType.__members__.values():
      raise ValueError(f"Invalid compression type {value}")
    # Should check value
    self._compression_type = value

  def __del__(self):
    self.close()

  def next(self) -> Optional[str]:
    """Read the next item from self._stream
    """

    nbytes = self._stream.read(4)
#   print(f"Read {len(nbytes)} bytes")
    nbytes = int.from_bytes(nbytes, "little", signed=False)
    if nbytes == 0:
      return None
#   print(f"Reading {nbytes}")
    data = self._stream.read(nbytes)
    if self._compression_type == CompressionType.ctype_uncompressed:
      return data

    if self._compression_type == CompressionType.ctype_snappy:
      return snappy.uncompress(data)

  def close(self):
    self._stream.close()
