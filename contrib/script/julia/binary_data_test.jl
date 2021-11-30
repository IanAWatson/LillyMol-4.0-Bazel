using Test
using BinaryData
using StringViews

function one_item()
  fname = tempname()
  writer::BinaryDataFileWriter = BinaryDataFileWriter()
  open(writer, fname) || return false
  written = "hello world"
  write(writer, written) || return false
  close(writer)

  reader::BinaryDataFileReader = BinaryDataFileReader()
  open(reader, fname) || return false
  returned = StringView(read(reader))
  return returned == written
end

function several_items()
  nitems = 100
  items = ["Item $(i) xyz" for i in 1:nitems]
  fname = tempname()
  writer = BinaryDataFileWriter()
  open(writer, fname) || return false
  for item in items
    write(writer, item) || return false
  end
  close(writer)

  reader = BinaryDataFileReader()
  open(reader, fname) || return false

  for i in 1:nitems
    returned = read(reader)
    StringView(returned) == items[i] || return false
  end

  return read(reader) == nothing
end

function snappy_compression()::Bool
  nitems = 100
  items = ["Item $(i) xyz" ^ 1000 for i in 1:nitems]
  fname = tempname()
  writer = BinaryDataFileWriter()
  set_compression(writer, BinaryData.snappy)
  open(writer, fname) || return false
  for item in items
    write(writer, item) || return false
  end
  close(writer)

  reader = BinaryDataFileReader()
  open(reader, fname) || return false

  for i in 1:nitems
    returned = read(reader)
    StringView(returned) == items[i] || return false
  end
  return read(reader) == nothing
end

function reader_cannot_open()
  reader = BinaryDataFileReader()
  return ! open(reader, "/this does not exist")
end

function writer_cannot_open()
  reader = BinaryDataFileWriter()
  return ! open(reader, "/this does not exist")
end

@test reader_cannot_open()
@test writer_cannot_open()
@test one_item()
@test several_items()
@test snappy_compression()
