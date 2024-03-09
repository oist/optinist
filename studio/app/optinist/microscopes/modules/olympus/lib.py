"""Olympus IDA wrapper module

* Utilities for IDA API

"""
import ctypes as ct

from studio.app.optinist.microscopes.modules.olympus.h_ida import (
    IDA_AXIS_INFO,
    IDA_PARAM,
    IDA_PARAM_ELEMENT,
    IDA_VALUE,
    IDA_AxisType,
    IDA_Result,
)

ida = None


def load_library(library_path: str):
    global ida
    ida = ct.cdll.LoadLibrary(library_path)
    return ida


def get_property_value(hAccessor, hProp, propName):
    size = ct.c_int()
    result = ida.GetPropertyValue(hAccessor, hProp, propName, None, 0, ct.byref(size))
    ret = None
    if result == IDA_Result.IDA_RESULT_SUCCESS:
        ret = (IDA_VALUE * size.value)()
        ida.GetPropertyValue(hAccessor, hProp, propName, ret, size, ct.byref(size))
    else:
        raise Exception(result)
    return result, ret


def get_frame_property(hAccessor, hImage, key, axisName=None, axisValue=None):
    hProp = ct.c_void_p()
    Param_AxisPosition = IDA_PARAM()
    if axisName:
        c_axisValue = ct.c_wchar_p(axisValue)
        Param_AxisPosition.nSize = 1
        element = IDA_PARAM_ELEMENT()
        element.pszKey = axisName
        element.pValue = ct.cast(c_axisValue, ct.c_void_p)
        Param_AxisPosition.pElements.contents = element
    result = ida.GetFrameProperty(
        hAccessor, hImage, key, ct.byref(Param_AxisPosition), ct.byref(hProp)
    )
    return result, hProp


def get_area_property(hAccessor, hArea, key, element_data=None):
    param = IDA_PARAM()
    if element_data:
        element = IDA_PARAM_ELEMENT()
        element.pszKey = element_data[0]
        element.pValue = element_data[1]

        param.pElements.contents = element
        param.nSize = 1
    else:
        param.nSize = 0
    hProp = ct.c_void_p()
    result = ida.GetAreaProperty(
        hAccessor, hArea, key, ct.byref(param), ct.byref(hProp)
    )
    return result, hProp


def get_lut(hAccessor, hArea, channel_id):
    size = ct.c_int()
    ida.GetLUT(hAccessor, hArea, channel_id, None, None, None, 0, ct.byref(size))
    pLUTR = (ct.c_int * size.value)()
    pLUTG = (ct.c_int * size.value)()
    pLUTB = (ct.c_int * size.value)()
    result = ida.GetLUT(
        hAccessor, hArea, channel_id, pLUTR, pLUTG, pLUTB, size, ct.byref(size)
    )
    return result, pLUTR, pLUTG, pLUTB


def get_image_body(hAccessor, hImage, rect):
    image_size = ct.c_uint64()
    image_buffer = None
    result = ida.GetImageBody(hAccessor, hImage, rect, None, None, ct.byref(image_size))
    if result != IDA_Result.IDA_RESULT_SUCCESS:
        raise Exception("Error: GetImageBody")

    if image_size.value != 0:
        image_buffer = (ct.c_uint8 * image_size.value)()
        result = ida.GetImageBody(
            hAccessor, hImage, rect, image_buffer, image_size, ct.byref(image_size)
        )
    return image_buffer


def get_image_axis(hAccessor, hImage):
    num_of_frame_axis = ct.c_int()
    result = ida.GetImageAxis(
        hAccessor, hImage, None, None, ct.byref(num_of_frame_axis)
    )
    if result != IDA_Result.IDA_RESULT_SUCCESS:
        raise Exception("Error: GetImageAxis")

    pFrameAxes = (IDA_AXIS_INFO * num_of_frame_axis.value)()

    result = ida.GetImageAxis(
        hAccessor, hImage.pFrameAxes, num_of_frame_axis, ct.byref(num_of_frame_axis)
    )
    if result != IDA_Result.IDA_RESULT_SUCCESS:
        raise Exception("Error: GetImageAxis")

    return pFrameAxes


def set_frame_axis_index(
    nLIndex, nZIndex, nTIndex, pRoiCollection, pAxisInfo, pAxes, pnAxisCount
):
    """
    // *****************************************************************************
    //  @brief  SetFrameAxisIndex
    //  @param[in]  nLIndex -- LAMBDA axis index if LAMBDA axis doesn't exist, set 0
    //  @param[in]  nZIndex -- ZSTACK axis index if ZSTACK axis doesn't exist, set 0
    //  @param[in]  nTIndex -- TIMELASE axis index if TIMELASE axis doesn't exist, set 0
    //  @param[in]	pAxisInfo- Axis Information getting from CAxisInfo class
    //  @param[out]  pAes
    //      -- This output is used by GetImage Function( CFrameManager::CFrameManager )
    //  @retval		void
    //  @note
    // *****************************************************************************
    """

    # 0: LAxis
    # 1: ZAxis
    # 2: TAxis
    KEY = ["LAMBDA", "ZSTACK", "TIMELAPSE"]
    nSize = [0] * 3
    bHasAxis = [pAxisInfo.exist(key) for key in KEY]
    pAxis = [pAxisInfo.get_axis(key) for key in KEY]
    pnAxisCount = 0
    for i in range(3):
        if bHasAxis[i]:
            nSize[i] = pAxis[i].get_max()
    if pRoiCollection.has_point_roi():  # Point
        pAxes[0].nNumber = nTIndex
        pAxes[0].nType = IDA_AxisType.IDA_AT_TIME
        pnAxisCount = 1
    elif pRoiCollection.has_multi_point_roi():  # Multipoint
        pAxes[0].nNumber = nTIndex
        pAxes[0].nType = IDA_AxisType.IDA_AT_TIME
        pnAxisCount = 1
    elif pRoiCollection.has_mapping_roi():  # Mapping
        if not bHasAxis[1]:
            pAxes[0].nNumber = nTIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 1
        else:
            pAxes[0].nNumber = nZIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_Z
            pAxes[1].nNumber = nTIndex
            pAxes[1].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 2
    elif pRoiCollection.has_line_roi():  # Line
        if not bHasAxis[0] and not bHasAxis[1] and bHasAxis[2]:  # XT
            pAxes[0].nNumber = nTIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 1
        elif not bHasAxis[0] and bHasAxis[1] and not bHasAxis[2]:  # XZ
            pAxes[0].nNumber = nZIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_Z
            pnAxisCount = 1
        elif not bHasAxis[0] and bHasAxis[1] and bHasAxis[2]:  # XZT
            pAxes[0].nNumber = nZIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_Z
            pAxes[1].nNumber = nTIndex
            pAxes[1].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 2
        else:
            # Note: Not implemented.
            pass
    else:  # XY
        if nSize[0] != 0 and nSize[1] != 0 and nSize[2] != 0:  # XYLZT
            pAxes[0].nNumber = nLIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_LAMBDA
            pAxes[1].nNumber = nZIndex
            pAxes[1].nType = IDA_AxisType.IDA_AT_Z
            pAxes[2].nNumber = nTIndex
            pAxes[2].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 3
        elif nSize[0] != 0 and nSize[1] != 0 and nSize[2] == 0:  # XYLZ
            pAxes[0].nNumber = nLIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_LAMBDA
            pAxes[1].nNumber = nZIndex
            pAxes[1].nType = IDA_AxisType.IDA_AT_Z
            pnAxisCount = 2
        elif nSize[0] != 0 and nSize[1] == 0 and nSize[2] != 0:  # XYLT
            pAxes[0].nNumber = nLIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_LAMBDA
            pAxes[1].nNumber = nTIndex
            pAxes[1].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 2
        elif nSize[0] == 0 and nSize[1] != 0 and nSize[2] != 0:  # XYZT
            pAxes[0].nNumber = nZIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_Z
            pAxes[1].nNumber = nTIndex
            pAxes[1].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 2
        elif nSize[0] != 0 and nSize[1] == 0 and nSize[2] == 0:  # XYL
            pAxes[0].nNumber = nLIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_LAMBDA
            pnAxisCount = 1
        elif nSize[0] == 0 and nSize[1] != 0 and nSize[2] == 0:  # XYZ
            pAxes[0].nNumber = nZIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_Z
            pnAxisCount = 1
        elif nSize[0] == 0 and nSize[1] == 0 and nSize[2] != 0:  # XYT
            pAxes[0].nNumber = nTIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 1
        elif nSize[0] == 0 and nSize[1] == 0 and nSize[2] == 0:  # XY
            pAxes[0].nNumber = nTIndex
            pAxes[0].nType = IDA_AxisType.IDA_AT_TIME
            pnAxisCount = 1
    del pAxis
    return pnAxisCount
