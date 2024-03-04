"""Olympus IDA wrapper module

* This module is not provided in IDA_Sample.

"""

import studio.app.optinist.microscopes.modules.olympus.lib as lib


class ScannerSettings:
    def __init__(self, hAccessor, hArea):
        self.__speed: float = None
        self.__roundtrip: str = None
        self.__pixelSpeed: float = None
        self.__lineSpeed: float = None
        self.__frameSpeed: float = None
        self.__seriesInterval: float = None
        self.__focus: int = None

        result, hProp = lib.get_area_property(hAccessor, hArea, "ScannerSettings")

        result, speed = lib.get_property_value(hAccessor, hProp, "speed")
        self.__speed = speed[0].value.dDouble
        del speed

        result, roundtrip = lib.get_property_value(hAccessor, hProp, "roundtrip")
        self.__roundtrip = roundtrip[0].value.pszString
        del roundtrip

        result, pixelSpeed = lib.get_property_value(hAccessor, hProp, "pixelSpeed")
        self.__pixelSpeed = pixelSpeed[0].value.dDouble
        del pixelSpeed

        result, lineSpeed = lib.get_property_value(hAccessor, hProp, "lineSpeed")
        self.__lineSpeed = lineSpeed[0].value.dDouble
        del lineSpeed

        result, frameSpeed = lib.get_property_value(hAccessor, hProp, "frameSpeed")
        self.__frameSpeed = frameSpeed[0].value.dDouble
        del frameSpeed

        result, seriesInterval = lib.get_property_value(
            hAccessor, hProp, "seriesInterval"
        )
        self.__seriesInterval = seriesInterval[0].value.dDouble
        del seriesInterval

        result, focus = lib.get_property_value(hAccessor, hProp, "focus")
        self.__focus = focus[0].value.nInteger
        del focus

        if hProp:
            lib.ida.ReleaseProperty(hAccessor, hProp)

    def print(self):
        print("Scanner Settings")
        print(f"\tspeed         = {self.__speed}")
        print(f"\troundtrip     = {self.__roundtrip}")
        print(f"\tpixelSpeed    = {self.__pixelSpeed}")
        print(f"\tlineSpeed     = {self.__lineSpeed}")
        print(f"\tframeSpeed    = {self.__frameSpeed}")
        print(f"\tseriesInterval= {self.__seriesInterval}")
        print(f"\tfocus         = {self.__focus}")

    def get_values(self):
        return {
            "speed": self.__speed,
            "roundtrip": self.__roundtrip,
            "pixel_speed": self.__pixelSpeed,
            "line_speed": self.__lineSpeed,
            "frame_speed": self.__frameSpeed,
            "series_interval": self.__seriesInterval,
            "focus": self.__focus,
        }

    @property
    def speed(self):
        return self.__speed

    @property
    def roundtrip(self):
        return self.__roundtrip

    @property
    def pixel_speed(self):
        return self.__pixelSpeed

    @property
    def line_speed(self):
        return self.__lineSpeed

    @property
    def frame_speed(self):
        return self.__frameSpeed

    @property
    def series_interval(self):
        return self.__seriesInterval

    @property
    def focus(self):
        return self.__focus
