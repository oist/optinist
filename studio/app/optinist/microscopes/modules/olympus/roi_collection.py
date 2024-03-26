"""Olympus IDA wrapper module

* Porting of IDA_Sample/RoiCollection.h,cpp

"""
import ctypes as ct

import studio.app.optinist.microscopes.modules.olympus.h_ida as h_ida
import studio.app.optinist.microscopes.modules.olympus.lib as lib


class Roi:
    def __init__(self):
        self.m_szId = None
        self.m_szName = None
        self.m_szType = None
        self.m_szShape = None

        self.m_dPanX = 0.0
        self.m_dPanY = 0.0
        self.m_dZoom = 0.0
        self.m_dZPosition = 0.0
        self.m_dRotation = 0.0
        self.m_vecPoints = []

    def set_id(self, pSrc):
        self.m_szId = pSrc
        return True

    def set_name(self, pSrc):
        self.m_szName = pSrc
        return True

    def set_type(self, pSrc):
        self.m_szType = pSrc
        return True

    def st_shape(self, pSrc):
        self.m_szShape = pSrc
        return True

    def set_rotation(self, pSrc):
        self.m_dRotation = pSrc
        return True

    def set_points(self, pSrc, nNumOfSrc):
        for p in pSrc:
            buf = h_ida.IDA_VALUE.IDA_POINT()
            buf.x = p.value.point.x
            buf.y = p.value.point.y
            self.m_vecPoints.append(buf)

        return True

    def set_pan_x(self, dPanX):
        self.m_dPanX = dPanX

    def set_pan_y(self, dPanY):
        self.m_dPanY = dPanY

    def set_zoom(self, dZoom):
        self.m_dZoom = dZoom

    def set_z(self, dZ):
        self.m_dZPosition = dZ

    def get_id(self):
        return self.m_szId

    def get_name(self):
        return self.m_szName

    def get_type(self):
        return self.m_szType

    def get_shape(self):
        return self.m_szShape

    def get_rotation(self):
        return self.m_dRotation

    def get_points(self):
        return self.m_vecPoints

    def get_pan_x(self):
        return self.m_dPanX

    def get_pan_y(self):
        return self.m_dPanY

    def get_zoom(self):
        return self.m_dZoom

    def get_z(self):
        return self.dZ

    def write_roi(self):
        print(f"\tId = {self.get_id()}")
        print(f"\t\tName = {self.get_name()}")
        print(f"\t\tShape = {self.get_shape()}")
        print(f"\t\tType = {self.get_type()}")
        print(f"\t\tRotation = {self.get_rotation()}")

        for p in self.get_points():
            print(f"\t\t\tx={p.x}, y={p.y}")

        print(f"\t\tPanX = {self.get_pan_x()}")
        print(f"\t\tPanY = {self.get_pan_y()}")
        print(f"\t\tZoom = {self.get_zoom()}")
        print(f"\t\tZ = {self.get_z()}")


class RoiCollection:
    def __init__(self, hAccessor, hArea, pKey_AnalysisROIList, pKey_AnalysisROIInfo):
        self.m_vecRois = []
        result, hPropList = lib.get_area_property(
            hAccessor, hArea, pKey_AnalysisROIList
        )
        result, pAnalysisROIDs = lib.get_property_value(hAccessor, hPropList, "id")
        for a_roi in pAnalysisROIDs:
            roi = Roi()

            # Get Image ROI Info from ID
            result, hPropInfo = lib.get_area_property(
                hAccessor,
                hArea,
                pKey_AnalysisROIInfo,
                ["roiId", ct.c_wchar_p(a_roi.value.pszString)],
            )
            # Get Analysis ROI Name
            result, pAnalysisROIName = lib.get_property_value(
                hAccessor, hPropInfo, "name"
            )
            roi.set_id(ct.p.value.pszString)
            roi.set_name(pAnalysisROIName[0].value.pszString)
            del pAnalysisROIName

            # Get Analysis ROI Type
            result, pAnalysisROIType = lib.get_property_value(
                hAccessor, hPropInfo, "type"
            )
            roi.set_type(pAnalysisROIType[0].value.pszString)
            del pAnalysisROIType

            # Get Analysis ROI Shape
            result, pAnalysisROIShape = lib.get_property_value(
                hAccessor, hPropInfo, "shape"
            )
            roi.set_shape(pAnalysisROIShape[0].value.pszString)
            del pAnalysisROIShape

            # Get Analysis ROI Rotation
            result, pAnalysisROIRotation = lib.get_property_value(
                hAccessor, hPropInfo, "rotation"
            )
            roi.set_rotation(pAnalysisROIRotation[0].value.dDouble)
            del pAnalysisROIRotation

            # Get Analysis ROI Data
            result, pAnalysisROIData = lib.get_property_value(
                hAccessor, hPropInfo, "data"
            )
            roi.set_points(pAnalysisROIData, -1)
            del pAnalysisROIData

            # Multipoint Data
            if roi.get_type() == "MULTI_POINT":
                # PanX
                result, pPanX = lib.get_property_value(hAccessor, hPropInfo, "panX")
                roi.set_pan_x(pPanX[0].value.dDouble)
                del pPanX

                # PanY
                result, pPanY = lib.get_property_value(hAccessor, hPropInfo, "panY")
                roi.set_pan_y(pPanY[0].value.dDouble)
                del pPanY

                # Zoom
                result, pZoom = lib.get_proiperty_value(hAccessor, hPropInfo, "zoom")
                roi.set_zoom(pZoom[0].value.dDouble)
                del pZoom

                # Z
                result, pZ = lib.get_property_value(hAccessor, hPropInfo, "zPosition")
                roi.set_z(pZ[0].value.dDouble)
                del pZ

            self.m_vecRois.append(roi)
            if hPropInfo:
                lib.ida.ReleaseProperty(hAccessor, hPropInfo)
        del pAnalysisROIDs
        if hPropList:
            lib.ida.ReleaseProperty(hAccessor, hPropList)

    def print(self):
        for r in self.m_vecRois:
            r.write_roi()

    def has_line_roi(self):
        for p in self.m_vecRois:
            if p.get_type() == "LineROI" or p.get_type() == "FreeLineROI":
                return True
        return False

    def has_point_roi(self):
        for r in self.m_vecRois:
            if r.get_type() == "PointROI":
                return True
        return False

    def has_multi_point_roi(self):
        for r in self.m_vecRois:
            if r.get_type() == "MultiPointROI":
                return True
        return False

    def has_mapping_roi(self):
        for r in self.m_vecRois:
            if r.get_type() == "MappingROI":
                return True
        return False
