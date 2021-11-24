mutable struct BinaryDataFileWriter
  stream::Union{IO, Nothing}
end  # struct BinaryDataFileWriter

BinaryDataFileWriter() = BinaryDataFileWriter(nothing)

MAGIC_NUMBER = 3177520567

function write_number(writer::BinaryDataFileWriter, value::Int)::Bool
  number::UInt32 = value
  number = hton(number)
# println("WRiting $(number) from $(value)")
  return write(writer.stream, number) == 4
end

function open(writer::BinaryDataFileWriter, fname::String)::Bool
  try
    writer.stream = open(fname, "w")
  catch SystemError
    return false
  end
  return write_number(writer, MAGIC_NUMBER)
end

function flush(writer::BinaryDataFileWriter)
  flush(writer.stream)
end

function close(writer::BinaryDataFileWriter)
  writer.stream == nothing && return
  close(writer.stream)
end

function write(writer::BinaryDataFileWriter, data::Any)::Bool
  write_number(writer, length(data))
  return write(writer.stream, data) == length(data)
end

mutable struct BinaryDataFileReader
  stream::Union{IO, Nothing}
end  # struct BinaryDataFileReader

BinaryDataFileReader() = BinaryDataFileReader(nothing)

function open(reader::BinaryDataFileReader, fname::String)::Bool
  try
    reader.stream = open(fname, "r")
  catch SystemError
    return false
  end
# println("OPened $(reader.stream):")
  magic = read_number(reader)
  return magic == MAGIC_NUMBER
end

function read_number(reader::BinaryDataFileReader)::UInt32
  try
    result = read(reader.stream, UInt32)
    return ntoh(result)
  catch EOFError
    return typemax(UInt32)
  end
end

function read(reader::BinaryDataFileReader)::Union{Nothing, Any}
  bytes = read_number(reader)
  bytes == typemax(UInt32) && return nothing
# println("Reading $(bytes) bytes")
  bytes == 0 && return nothing
  return read(reader.stream, bytes)
end

function close(reader::BinaryDataFileReader)
  reader.stream == nothing && return
  close(reader.stream)
end
