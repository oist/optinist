"""Olympus IDA wrapper module

* Porting of IDA_Sample/AxisInfo.h,cpp

"""
import ctypes as ct

import studio.app.optinist.microscopes.modules.olympus.h_ida as h_ida
import studio.app.optinist.microscopes.modules.olympus.lib as lib


class AxisIndex:
    def __init__(self):
        self.m_bIsExit = False
        self.m_nIndex = 0
        self.m_nType = h_ida.IDA_AxisType.IDA_AT_TIME

    def set_exit(self, bIsExit):
        self.m_bIsExit = bIsExit

    def set_index(self, nIndex):
        self.m_nIndex = nIndex

    def set_type(self, nType):
        self.m_nType = nType

    def get_exit(self):
        return self.m_bIsExit

    def get_index(self):
        return self.m_nIndex

    def get_type(self):
        return self.m_nType


class AxisPosition:
    def __init__(self):
        self.m_bIsExit = False
        self.m_dPos = 0.0
        self.m_nType = h_ida.IDA_AxisType.IDA_AT_TIME

    def set_exit(self, bIsExit):
        self.m_bIsExit = bIsExit

    def set_position(self, dPos):
        self.m_dPos = dPos

    def set_type(self, nType):
        self.m_nType = nType

    def get_exit(self):
        return self.m_bIsExit

    def get_position(self):
        return self.m_dPos

    def get_type(self):
        return self.m_nType


class Axis:
    def __init__(self):
        self.m_dStart = 0.0
        self.m_dEnd = 0.0
        self.m_dStep = 0.0
        self.m_nMax = 0

    def set_start(self, val):
        self.m_dStart = val

    def set_end(self, val):
        self.m_dEnd = val

    def set_step(self, val):
        self.m_dStep = val

    def set_max(self, val):
        self.m_nMax = val

    def get_start(self):
        return self.m_dStart

    def get_end(self):
        return self.m_dEnd

    def get_step(self):
        return self.m_dStep

    def get_max(self):
        return self.m_nMax


class AxisInfo:
    def __init__(self, hAccessor, hArea):
        self.m_axes = {}

        result, hPropAxes = lib.get_area_property(hAccessor, hArea, "Axes")
        result, pAxes = lib.get_property_value(hAccessor, hPropAxes, "axisName")

        for c_axis in pAxes:
            axis = Axis()
            result, hPropAxis = lib.get_area_property(
                hAccessor,
                hArea,
                "AxisInfo",
                [
                    "axisName",
                    ct.cast(ct.c_wchar_p(c_axis.value.pszString), ct.c_void_p),
                ],
            )
            if result == h_ida.IDA_Result.IDA_RESULT_SUCCESS:
                result, pStart = lib.get_property_value(hAccessor, hPropAxis, "start")
                axis.set_start(pStart[0].value.dDouble)
                del pStart

                result, pEnd = lib.get_property_value(hAccessor, hPropAxis, "end")
                axis.set_end(pEnd[0].value.dDouble)
                del pEnd

                result, pStep = lib.get_property_value(hAccessor, hPropAxis, "step")
                axis.set_step(pStep[0].value.dDouble)
                del pStep

                result, pMax = lib.get_property_value(hAccessor, hPropAxis, "maxSize")
                axis.set_max(pMax[0].value.nInteger)
                del pMax
                self.m_axes[c_axis.value.pszString] = axis
                if hPropAxis:
                    lib.ida.ReleaseProperty(hAccessor, hPropAxis)
        if pAxes:
            del pAxes

        if hPropAxes:
            lib.ida.ReleaseProperty(hAccessor, hPropAxes)

    def get_axis(self, name):
        return self.m_axes.get(name, None)

    def exist(self, name):
        return name in self.m_axes

    def print(self):
        print("Axis Information")
        for name, axis in self.m_axes.items():
            print(f"\taxisName = {name}")
            print(f"\t\tstart = {axis.get_start()}")
            print(f"\t\tstep = {axis.get_step()}")
            print(f"\t\tend = {axis.get_end()}")
            print(f"\t\tmax = {axis.get_max()}")

    def get_values(self):
        result = {}
        for name, axis in self.m_axes.items():
            result[name] = {
                "start": axis.get_start(),
                "step": axis.get_step(),
                "end": axis.get_end(),
                "max": axis.get_max(),
            }
        return result
