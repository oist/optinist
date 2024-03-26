"""Olympus IDA wrapper module

* Porting of IDA_Sample/SystemInfo.h,cpp

"""
import studio.app.optinist.microscopes.modules.olympus.lib as lib


class SystemInfo:
    def __init__(self, hAccessor, hArea):
        # variables
        self.m_szSystemName = None
        self.m_szSystemVersion = None
        self.m_szDeviceName = None
        self.m_szUserName = None

        result, hProp = lib.get_area_property(hAccessor, hArea, "SystemInfo")
        result, pSysName = lib.get_property_value(hAccessor, hProp, "name")
        self.m_szSystemName = pSysName[0].value.pszString
        del pSysName

        result, pSysVer = lib.get_property_value(hAccessor, hProp, "version")
        self.m_szSystemVersion = pSysVer[0].value.pszString
        del pSysVer

        result, pDeviceName = lib.get_property_value(hAccessor, hProp, "deviceName")
        self.m_szDeviceName = pDeviceName[0].value.pszString
        del pDeviceName

        result, pUserName = lib.get_property_value(hAccessor, hProp, "userName")
        self.m_szUserName = pUserName[0].value.pszString
        del pUserName

        if hProp:
            lib.ida.ReleaseProperty(hAccessor, hProp)

    def print(self):
        print("System Information")
        print(f"\tname = {self.m_szSystemName}")
        print(f"\tversion = {self.m_szSystemVersion}")
        print(f"\tdeviceName = {self.m_szDeviceName}")
        print(f"\tuserName = {self.m_szUserName}")

    def get_values(self):
        return {
            "system_name": self.m_szSystemName,
            "system_version": self.m_szSystemVersion,
            "device_name": self.m_szDeviceName,
            "user_name": self.m_szUserName,
        }

    @property
    def system_name(self):
        return self.m_szSystemName

    @property
    def system_version(self):
        return self.m_szSystemVersion

    @property
    def device_name(self):
        return self.m_szDeviceName

    @property
    def user_name(self):
        return self.m_szUserName
