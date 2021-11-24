# Tester for binary_data_file

from absl.testing import absltest
from absl.testing import parameterized

import binary_file

class TestBinary(absltest.TestCase):
  def testOneItem(self):
    fname = self.create_tempfile()
    writer = binary_file.BinaryDataFileWriter()
    self.assertTrue(writer.open(fname))
    data = "hello world"
    self.assertTrue(writer.write(data))
    writer.close()

    reader = binary_file.BinaryDataFileReader(fname)
    returned = reader.next()
    self.assertIsNotNone(returned)
    self.assertEqual(data, returned.decode("utf-8"))

    self.assertIsNone(reader.next())
    

  def testMultiple(self):
    """Test many items"""
    nitems = 100
    items = [f"foo{i}" for i in range(0, nitems)]
    fname = self.create_tempfile()
    writer = binary_file.BinaryDataFileWriter()
    self.assertTrue(writer.open(fname))
    for item in items:
      self.assertTrue(writer.write(item))
    writer.close()

    reader = binary_file.BinaryDataFileReader(fname)
    for i in range(0, nitems):
      returned = reader.next()
      self.assertEqual(items[i], returned.decode("utf-8"))

    self.assertIsNone(reader.next())

if __name__ == "__main__":
  absltest.main()
