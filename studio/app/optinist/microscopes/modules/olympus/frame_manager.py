import ctypes as ct

import studio.app.optinist.microscopes.modules.olympus.h_ida as h_ida
import studio.app.optinist.microscopes.modules.olympus.lib as lib
from studio.app.optinist.microscopes.modules.olympus.axis_info import (
    AxisIndex,
    AxisPosition,
)
from studio.app.optinist.microscopes.modules.olympus.roi_collection import Roi


class FrameManager:
    def __init__(self, hAccessor, hArea, pszChannelId, pAxes, nNumOfAxes):
        self.m_hAccessor = hAccessor
        self.m_hArea = hArea
        self.m_pszChannelId = pszChannelId
        self.m_nAxisCount = nNumOfAxes
        self.m_pAxes = []
        for axis in pAxes:
            ainfo = h_ida.IDA_AXIS_INFO()
            ainfo.nNumber = axis.nNumber
            ainfo.nType = axis.nType
            self.m_pAxes.append(ainfo)
        self.m_hImage = ct.c_void_p()
        self.m_rect = None
        self.m_pucImageBuffer = None
        self.m_pucImageBuffer_asWORD = None
        self.m_vecAxisIndex = []
        self.m_vecAxisPosition = []
        self.m_vecRois = []
        lib.ida.GetImage(
            hAccessor, hArea, pszChannelId, pAxes, nNumOfAxes, ct.byref(self.m_hImage)
        )

    def get_image_body(self, rect):
        self.m_rect = rect
        self.m_pucImageBuffer = lib.get_image_body(
            self.m_hAccessor, self.m_hImage, ct.byref(self.m_rect)
        )
        return self.m_pucImageBuffer

    def release_image_body(self):
        if self.m_pucImageBuffer:
            self.m_pucImagebuffer = None
            self.m_pucImageBuffer_asWORD = ct.cast(
                self.m_pucImageBuffer, ct.POINTER(ct.c_uint16)
            )

    def write_image_body(self, filename):
        # TODO: Do something here
        pass

    def write_image_body_binary(self, filename):
        # TODO: Do something here
        pass

    def get_pixel_value_tm(self, myDataCnt):
        return self.m_pucImageBuffer_asWORD[myDataCnt]

    def pucBuffer_to_WORD_TM(self, width, height):
        # TODO: The width and height should be obtained externally as appropriate.
        #   (If it is difficult to obtain them here, they should be obtained elsewhere.)
        buffer_size = ct.c_uint16 * width * height

        self.m_pucImageBuffer_asWORD = buffer_size.from_buffer(self.m_pucImageBuffer)

        return self.m_pucImageBuffer_asWORD

    def get_frame_index(self):
        pFrameAxes = lib.get_image_axis(self.m_hAccessor, self.m_hImage)
        for f in pFrameAxes:
            axis_index = AxisIndex()
            axis_index.set_exit(True)
            # TODO: Is the following code correct? (p.nType, p.nNumber)
            axis_index.set_type(ct.p.nType)
            axis_index.set_index(ct.p.nNumber)
            self.m_vecAxisIndex.append(axis_index)

        del pFrameAxes
        return self.m_vecAxisIndex

    def write_frame_index(self):
        for ai in self.m_vecAxisIndex:
            if ai.get_exit():
                print(f"\tType={ai.get_type()}, Number={ai.get_index()}")

    def get_frame_position(self):
        result, hProp = lib.get_frame_property(
            self.m_hAccessor, self.m_hImage, "AxisPosition", "axisName", "TIMELAPSE"
        )
        if result == 0:
            result, pAxisPosition = lib.get_property_value(
                self.m_hAccessor, hProp, "position"
            )
            axis_pos = AxisPosition()
            axis_pos.set_type(h_ida.IDA_AxisType.IDA_AT_TIME)
            axis_pos.set_exit(True)
            axis_pos.set_position(pAxisPosition[0].value.dDouble)
            self.m_vecAxisPosition.append(axis_pos)
            del pAxisPosition
        if hProp:
            lib.ida.ReleaseProperty(self.m_hAccessor, hProp)

        result, hProp = lib.get_frame_property(
            self.m_hAccessor, self.m_hImage, "AxisPosition", "axisName", "ZSTACK"
        )
        if result == 0:
            result, pAxisPosition = lib.get_property_value(
                self.m_hAccessor, hProp, "position"
            )
            axis_pos = AxisPosition()
            axis_pos.set_type(h_ida.IDA_AxisType.IDA_AT_Z)
            axis_pos.set_exit(True)
            axis_pos.set_position(pAxisPosition[0].value.dDouble)
            self.m_vecAxisPosition.append(axis_pos)
            del pAxisPosition
        if hProp:
            lib.ida.ReleaseProperty(self.m_hAccessor, hProp)

        result, hProp = lib.get_frame_property(
            self.m_hAccessor, self.m_hImage, "AxisPosition", "axisName", "LAMBDA"
        )
        if result == 0:
            result, pAxisPosition = lib.get_property_value(
                self.m_hAccessor, hProp, "position"
            )
            axis_pos = AxisPosition()
            axis_pos.set_type(h_ida.IDA_AxisType.IDA_AT_LAMBDA)
            axis_pos.set_exit(True)
            axis_pos.set_position(pAxisPosition[0].value.dDouble)
            self.m_vecAxisPosition.append(axis_pos)
            del pAxisPosition
        if hProp:
            lib.ida.ReleaseProperty(hProp)

        return self.m_vecAxisPosition

    def write_frame_position(self):
        for ap in self.m_vecAxisPosition:
            if ap.get_exit():
                if ap.get_type() == h_ida.IDA_AxisType.IDA_AT_LAMBDA:
                    print(f"\tLAMBDA={ap.get_position()}")
                elif ap.get_type() == h_ida.IDA_AxisType.IDA_AT_Z:
                    print(f"\tZSTACK={ap.get_position()}")
                elif ap.get_type() == h_ida.IDA_AxisType.IDA_AT_TIME:
                    print(f"\tTIMELAPSE={ap.get_position()}")

    def get_timestamp_channel_tm(self):
        my_timestamp_channel = 0
        for idx, ap in enumerate(self.m_vecAxisPosition):
            if ap.get_exit():
                if ap.get_type() == h_ida.IDA_AxisType.IDA_AT_TIME:
                    my_timestamp_channel = idx
        return my_timestamp_channel

    def get_timestamp_tm(self, mychannel):
        mytimestamp = self.m_vecAxisPosition[mychannel].get_position()
        return mytimestamp

    def get_frame_roi(self):
        result, hProp = lib.get_frame_property(
            self.m_hAccessor, self.m_hImage, "StimulationROIList"
        )
        result, pAnalysisROIDs = lib.get_property_value(self.m_hAccessor, hProp, "id")
        for aroi in pAnalysisROIDs:
            roi = Roi()
            # Get Image ROI Info from ID
            result, hPropInfo = lib.get_frame_property(
                self.m_hAccessor,
                self.m_hImage,
                "StimulationROIInfo",
                "roiId",
                ct.c_wchar_p(aroi.value.pszString),
            )
            # Get Analysis ROI Name
            result, pAnalysisROIName = lib.get_property_value(
                self.m_hAccessor, hPropInfo, "name"
            )
            roi.set_id(aroi.value.pszString)
            roi.set_name(pAnalysisROIName[0].value.pszString)
            del pAnalysisROIName

            # Get Analysis ROI Type
            result, pAnalysisROIType = lib.get_property_value(
                self.m_hAccessor, hPropInfo, "type"
            )
            roi.set_type(pAnalysisROIType[0].value.pszString)
            del pAnalysisROIType

            # Get Analysis ROI Shape
            result, pAnalysisROIShape = lib.get_property_value(
                self.m_hAccessor, hPropInfo, "shape"
            )
            roi.set_shape(pAnalysisROIShape[0].value.pszString)
            del pAnalysisROIShape

            # Get Analysis ROI Rotation
            result, pAnalysisROIRotation = lib.get_property_value(
                self.m_hAccessor, hPropInfo, "rotation"
            )
            roi.set_rotation(pAnalysisROIRotation[0].value.dDouble)
            del pAnalysisROIRotation

            # Get Analysis ROI Data
            result, pAnalysisData = lib.get_property_value(
                self.m_hAccessor, hPropInfo, "data"
            )

            # TODO: linterで警告が生じているため、コメントアウト（別途修正）
            # roi.set_points(pAnalysisROIData, -1)

            del pAnalysisData

            if roi.get_type() == "MULTI_POINT":
                # PanX
                result, pPanX = lib.get_property_value(
                    self.m_hAccessor, hPropInfo, "panX"
                )
                roi.set_pan_x(pPanX[0].value.dDouble)
                del pPanX

                # PanY
                result, pPanY = lib.get_property_value(
                    self.m_hAccessor, hPropInfo, "panY"
                )
                roi.set_pan_y(pPanY[0].value.dDouble)
                del pPanY

                # Zoom
                result, pZoom = lib.get_property_value(
                    self.m_hAccessor, hPropInfo, "zoom"
                )
                roi.set_zoom(pZoom[0].value.dDouble)
                del pZoom

                # Z
                result, pZ = lib.get_property_value(
                    self.m_hAccessor, hPropInfo, "zPosition"
                )
                roi.set_z(pZ[0].value.dDouble)
                del pZ

            self.m_vecRois.append(roi)
            if hPropInfo:
                lib.ida.ReleaseProperty(self.m_hAccessor, hPropInfo)
        del pAnalysisROIDs
        if hProp:
            lib.ida.ReleaseProperty(self.m_hAccessor, hProp)

        return self.m_vecRois

    def write_frame_roi(self):
        for r in self.m_vecRois:
            r.write_roi()
