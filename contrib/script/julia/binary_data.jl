using Snappy
using StringViews

@enum CompressionType uncompressed=0 snappy=1 invalid=2

mutable struct BinaryDataFileWriter
  stream::Union{IO, Nothing}
  compression_type::CompressionType
end  # struct BinaryDataFileWriter

BinaryDataFileWriter() = BinaryDataFileWriter(nothing, uncompressed)

MAGIC_NUMBER = 3177520567

function set_compression(writer::BinaryDataFileWriter, ctype::CompressionType)::Bool
  if ! isnothing(writer.stream)
    println("BinaryDataFileWriter:set_compression:already open")
    return false
  end
  writer.compression_type = ctype
  true
end

function write_number(writer::BinaryDataFileWriter, value::Int)::Bool
  number::UInt32 = value
  number = hton(number)
# println("WRiting $(number) from $(value)")
  return write(writer.stream, number) == 4
end

function write_header(writer::BinaryDataFileWriter)::Bool
  write_number(writer, MAGIC_NUMBER) || return false
  write_number(writer, Int(writer.compression_type))
end

function open(writer::BinaryDataFileWriter, fname::String)::Bool
  try
    writer.stream = open(fname, "w")
  catch SystemError
    return false
  end
  return write_header(writer)
end

function flush(writer::BinaryDataFileWriter)
  flush(writer.stream)
end

function close(writer::BinaryDataFileWriter)
  writer.stream == nothing && return
  close(writer.stream)
end

function write(writer::BinaryDataFileWriter, data::Any)::Bool
  if writer.compression_type == uncompressed
    write_number(writer, length(data))
    return write(writer.stream, data) == length(data)
  end

  if writer.compression_type == snappy
    compressed = Snappy.compress(Vector{UInt8}(data))  # Likely make a copy...
    write_number(writer, length(compressed))
#   println("Data compressed to $(length(compressed)) bytes")
    return write(writer.stream, compressed) == length(compressed)
  end

  return false
end

mutable struct BinaryDataFileReader
  stream::Union{IO, Nothing}
  compression_type::CompressionType
end  # struct BinaryDataFileReader

BinaryDataFileReader() = BinaryDataFileReader(nothing, uncompressed)

function open(reader::BinaryDataFileReader, fname::String)::Bool
  try
    reader.stream = open(fname, "r")
  catch SystemError
    return false
  end
# println("OPened $(reader.stream):")
  return read_header(reader)
end

function read_header(reader::BinaryDataFileReader)::Bool
  MAGIC_NUMBER == read_number(reader) || return false
  number_to_compression = [uncompressed, snappy, invalid]
  reader.compression_type = number_to_compression[read_number(reader) + 1]
  return true
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
  data = read(reader.stream, bytes)
  if reader.compression_type == uncompressed
    return data
  end

  if reader.compression_type == snappy
#   println("uncompress $(typeof(data)) size $(length(data))")
    return Snappy.uncompress(data)
  end
end

function close(reader::BinaryDataFileReader)
  reader.stream == nothing && return
  close(reader.stream)
end
