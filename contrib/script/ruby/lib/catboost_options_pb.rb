# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: catboost_options.proto

require 'google/protobuf'

Google::Protobuf::DescriptorPool.generated_pool.build do
  add_file("catboost_options.proto", :syntax => :proto3) do
    add_message "CatboostOptions" do
      map :string_option, :string, :string, 1
      map :int_option, :string, :int32, 2
      map :float_option, :string, :float, 3
    end
  end
end

CatboostOptions = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("CatboostOptions").msgclass
