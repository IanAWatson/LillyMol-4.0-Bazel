# Interface to Foundational/data_source/binary_data_file.h

module BinaryData

import Base.finalizer

import Base.close
import Base.read
import Base.open
import Base.write

include("binary_data.jl")

export BinaryDataFileReader
export BinaryDataFileWriter
#export open(::BinaryDataFileWriter, ::String)
#export open
#export read
export set_compression
export CompressionType
export snappy

end # module BinaryData
