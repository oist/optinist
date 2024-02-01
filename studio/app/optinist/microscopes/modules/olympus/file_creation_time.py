"""Olympus IDA wrapper module

* Porting of IDA_Sample/FileCreationTime.h,cpp

"""
import studio.app.optinist.microscopes.modules.olympus.lib as lib


class FileCreationTime:
    def __init__(self, hAccessor, hArea):
        # Variables
        self.m_szCreationTime = None

        result, hProp = lib.get_area_property(hAccessor, hArea, "CreationDateTime")
        result, pCreationDateTime = lib.get_property_value(hAccessor, hProp, "dateTime")
        self.m_szCreationTime = pCreationDateTime[0].value.pszString
        del pCreationDateTime

        if hProp:
            lib.ida.ReleaseProperty(hAccessor, hProp)

    def print(self):
        print("File Creation Time")
        print(f"\tTime={self.m_szCreationTime}")

    def get_values(self):
        return {
            "creation_time": self.m_szCreationTime,
        }

    @property
    def creation_time(self):
        return self.m_szCreationTime
