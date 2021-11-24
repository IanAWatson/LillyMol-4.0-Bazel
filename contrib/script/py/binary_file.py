# Python tools to read and write files used by
# BinaryDataFileReader and BinaryDataFileWriter

import ctypes

from typing import Optional

MAGIC_NUMBER = 3177520567

class BinaryDataFileWriter:
  """Create files compatible with Foundational/datasource/binary_data_file.h
  """
  def __init__(self):
    self._stream = None

  def __del__(self):
    self.close()

  def open(self, fname: str) -> bool:
    """Open `fname` for writing.

    Args:
      fname: the file to open
    Returns:
      True on success
    """
    self._stream = open(fname, "wb")
    return self.write_number(MAGIC_NUMBER)

  def close(self):
    self._stream.close()

  def flush(self):
    self._stream.flush();

  def write(self, data):
    self.write_number(len(data))
    if len(data) > 0:
      self._stream.write(str.encode(data))
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
    self._stream = open(fname, "rb")
    value = self._stream.read(4)
    value = int.from_bytes(value, "little", signed=False)
    if value != MAGIC_NUMBER:
      raise ValueError(f"Invalid magic number got {value} not MAGIC_NUMBER")
    return

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
    return self._stream.read(nbytes)

  def close(self):
    self._stream.close()
