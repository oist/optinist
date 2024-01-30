"""Olympus IDA wrapper module

* Porting of IDA/include/Idaldll.h

"""
import ctypes as ct
from enum import IntEnum, auto


class IDA_Result(IntEnum):
    IDA_RESULT_SUCCESS = 0
    IDA_RESULT_INVALID_PARAM = auto()
    IDA_RESULT_ILLEGAL_STATE = auto()
    IDA_RESULT_ACCESSOR_NOT_FOUND = auto()
    IDA_RESULT_FILE_NOT_FOUND = auto()
    IDA_RESULT_PROPERTY_NOT_FOUND = auto()
    IDA_RESULT_CHANNEL_NOT_FOUND = auto()
    IDA_RESULT_AMBIGUOUS_FRAME = auto()
    IDA_RESULT_NO_LUT = auto()
    IDA_RESULT_UNSUPPORTED_OPERATION = auto()
    IDA_RESULT_NO_IMAGE = auto()
    IDA_RESULT_ALLOCATION_FAILED = auto()
    IDA_RESULT_FILE_CORRUPT = auto()
    IDA_RESULT_FAILURE = auto()
    IDA_INVALID_IMAGE = auto()
    IDA_RESULT_INSUFFICIENT_BUFFER = auto()
    IDA_RESULT_IMAGE_NOT_ACQUIRED = auto()
    IDA_RESULT_NO_MORE_IMAGE = auto()


class IDA_OpenMode(IntEnum):
    IDA_OM_READ = 0x00
    IDA_OM_WRITE = 0x01
    IDA_OM_NORMAL = 0x00
    IDA_OM_FORCE = 0x02
    IDA_OM_STREAMING = 0x04
    IDA_OM_READ_WRITE = 0x08
    IDA_OM_APPEND = 0x10


class CMN_RECT(ct.Structure):
    _fields_ = [
        ("x", ct.c_uint64),
        ("y", ct.c_uint64),
        ("width", ct.c_uint64),
        ("height", ct.c_uint64),
    ]


class IDA_POINT(ct.Structure):
    _fields_ = [
        ("x", ct.c_int),
        ("y", ct.c_int),
    ]


class IDA_VALUE_UN(ct.Union):
    _fields_ = [
        ("nInteger", ct.c_int),
        ("dDouble", ct.c_double),
        ("pszString", ct.c_wchar_p),
        ("point", IDA_POINT),
        ("rRect", CMN_RECT),
    ]


class IDA_VALUE(ct.Structure):
    _fields_ = [
        ("nType", ct.c_int),
        ("value", IDA_VALUE_UN),
    ]


class IDA_PARAM_ELEMENT(ct.Structure):
    _fields_ = [
        ("pszKey", ct.c_wchar_p),
        ("pValue", ct.c_void_p),
    ]


class IDA_PARAM(ct.Structure):
    _fields_ = [
        ("nSize", ct.c_int),
        ("pElements", ct.POINTER(IDA_PARAM_ELEMENT)),
    ]


class IDA_AXIS_INFO(ct.Structure):
    _fields_ = [
        ("nType", ct.c_int),
        ("nNumber", ct.c_int64),
    ]


class IDA_AxisType(IntEnum):
    IDA_AT_TIME = 0
    IDA_AT_Z = 1
    IDA_AT_LAMBDA = 2
