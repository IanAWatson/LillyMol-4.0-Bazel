# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: Utilities/General/class_label_translation.proto
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n/Utilities/General/class_label_translation.proto\x12\x15\x43lassLabelTranslation\"\xa0\x02\n\x15\x43lassLabelTranslation\x12O\n\nto_numeric\x18\x01 \x03(\x0b\x32;.ClassLabelTranslation.ClassLabelTranslation.ToNumericEntry\x12Q\n\x0b\x63lass_count\x18\x02 \x03(\x0b\x32<.ClassLabelTranslation.ClassLabelTranslation.ClassCountEntry\x1a\x30\n\x0eToNumericEntry\x12\x0b\n\x03key\x18\x01 \x01(\t\x12\r\n\x05value\x18\x02 \x01(\x05:\x02\x38\x01\x1a\x31\n\x0f\x43lassCountEntry\x12\x0b\n\x03key\x18\x01 \x01(\t\x12\r\n\x05value\x18\x02 \x01(\r:\x02\x38\x01\x62\x06proto3')



_CLASSLABELTRANSLATION = DESCRIPTOR.message_types_by_name['ClassLabelTranslation']
_CLASSLABELTRANSLATION_TONUMERICENTRY = _CLASSLABELTRANSLATION.nested_types_by_name['ToNumericEntry']
_CLASSLABELTRANSLATION_CLASSCOUNTENTRY = _CLASSLABELTRANSLATION.nested_types_by_name['ClassCountEntry']
ClassLabelTranslation = _reflection.GeneratedProtocolMessageType('ClassLabelTranslation', (_message.Message,), {

  'ToNumericEntry' : _reflection.GeneratedProtocolMessageType('ToNumericEntry', (_message.Message,), {
    'DESCRIPTOR' : _CLASSLABELTRANSLATION_TONUMERICENTRY,
    '__module__' : 'Utilities.General.class_label_translation_pb2'
    # @@protoc_insertion_point(class_scope:ClassLabelTranslation.ClassLabelTranslation.ToNumericEntry)
    })
  ,

  'ClassCountEntry' : _reflection.GeneratedProtocolMessageType('ClassCountEntry', (_message.Message,), {
    'DESCRIPTOR' : _CLASSLABELTRANSLATION_CLASSCOUNTENTRY,
    '__module__' : 'Utilities.General.class_label_translation_pb2'
    # @@protoc_insertion_point(class_scope:ClassLabelTranslation.ClassLabelTranslation.ClassCountEntry)
    })
  ,
  'DESCRIPTOR' : _CLASSLABELTRANSLATION,
  '__module__' : 'Utilities.General.class_label_translation_pb2'
  # @@protoc_insertion_point(class_scope:ClassLabelTranslation.ClassLabelTranslation)
  })
_sym_db.RegisterMessage(ClassLabelTranslation)
_sym_db.RegisterMessage(ClassLabelTranslation.ToNumericEntry)
_sym_db.RegisterMessage(ClassLabelTranslation.ClassCountEntry)

if _descriptor._USE_C_DESCRIPTORS == False:

  DESCRIPTOR._options = None
  _CLASSLABELTRANSLATION_TONUMERICENTRY._options = None
  _CLASSLABELTRANSLATION_TONUMERICENTRY._serialized_options = b'8\001'
  _CLASSLABELTRANSLATION_CLASSCOUNTENTRY._options = None
  _CLASSLABELTRANSLATION_CLASSCOUNTENTRY._serialized_options = b'8\001'
  _CLASSLABELTRANSLATION._serialized_start=75
  _CLASSLABELTRANSLATION._serialized_end=363
  _CLASSLABELTRANSLATION_TONUMERICENTRY._serialized_start=264
  _CLASSLABELTRANSLATION_TONUMERICENTRY._serialized_end=312
  _CLASSLABELTRANSLATION_CLASSCOUNTENTRY._serialized_start=314
  _CLASSLABELTRANSLATION_CLASSCOUNTENTRY._serialized_end=363
# @@protoc_insertion_point(module_scope)
