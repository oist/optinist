"""Olympus IDA wrapper module

* Porting of IDA_Sample/AreaImageSize.h,cpp

"""
import studio.app.optinist.microscopes.modules.olympus.lib as lib


class AreaImageSize:
    def __init__(self, hAccessor, hArea):
        result, hProp = lib.get_area_property(hAccessor, hArea, "ImageSize")
        result, pImageSize = lib.get_property_value(hAccessor, hProp, "size")
        self.m_nX = pImageSize[0].value.point.x
        self.m_nY = pImageSize[0].value.point.y
        if pImageSize:
            del pImageSize
        if hProp:
            lib.ida.ReleaseProperty(hAccessor, hProp)

    def get_x(self):
        return self.m_nX

    def get_y(self):
        return self.m_nY

    def print(self):
        print("Image Size")
        print(f"\tx={self.m_nX}")
        print(f"\ty={self.m_nY}")

    def get_values(self):
        return {
            "x": self.m_nX,
            "y": self.m_nY,
        }
