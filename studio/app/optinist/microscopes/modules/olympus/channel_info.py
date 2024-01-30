"""Olympus IDA wrapper module

* Porting of IDA_Sample/ChannelInfo.h,cpp

"""
import ctypes as ct

import studio.app.optinist.microscopes.modules.olympus.h_ida as h_ida
import studio.app.optinist.microscopes.modules.olympus.lib as lib


class ChannelInfo:
    def __init__(self, hAccessor, hArea):
        # Variables
        self.m_vecpszChannelIdList = []
        self.m_vecpszChannelName = []
        self.m_vecnChannelDepth = []
        self.m_vecnChannelBitCount = []
        self.m_vecpnChannelLUTR = []
        self.m_vecpnChannelLUTG = []
        self.m_vecpnChannelLUTB = []

        # Init
        result, hPropEnabled = lib.get_area_property(
            hAccessor, hArea, "EnabledChannelIdList"
        )
        result, pChIDs = lib.get_property_value(hAccessor, hPropEnabled, "id")
        # num_ch = len(pChIDs)

        # Get Channel Info
        for ch in pChIDs:
            self.m_vecpszChannelIdList.append(ch.value.pszString)

            result, hPropChannel = lib.get_area_property(
                hAccessor,
                hArea,
                "ChannelInfo",
                ["channelId", ct.cast(ct.c_wchar_p(ch.value.pszString), ct.c_void_p)],
            )

            if result == h_ida.IDA_Result.IDA_RESULT_SUCCESS:
                # Get Ch Names
                result, pChNames = lib.get_property_value(
                    hAccessor, hPropChannel, "name"
                )
                self.m_vecpszChannelName.append(pChNames[0].value.pszString)
                del pChNames

                # Get Ch depth (bit depth)
                result, pChDepth = lib.get_property_value(
                    hAccessor, hPropChannel, "depth"
                )
                self.m_vecnChannelDepth.append(pChDepth[0].value.nInteger)
                del pChDepth

                # Get Ch depth (available bit depth)
                result, pChAvaiDepth = lib.get_property_value(
                    hAccessor, hPropChannel, "bitCount"
                )
                self.m_vecnChannelBitCount.append(pChAvaiDepth[0].value.nInteger)
                del pChAvaiDepth

                # LUT
                result, lut_r, lut_g, lut_b = lib.get_lut(
                    hAccessor, hArea, ch.value.pszString
                )
                self.m_vecpnChannelLUTR.append(lut_r)
                self.m_vecpnChannelLUTG.append(lut_g)
                self.m_vecpnChannelLUTB.append(lut_b)

            if hPropChannel:
                lib.ida.ReleaseProperty(hAccessor, hPropChannel)

        if hPropEnabled:
            lib.ida.ReleaseProperty(hAccessor, hPropEnabled)
        if pChIDs:
            del pChIDs

    def get_num_of_channel(self):
        return len(self.m_vecpszChannelIdList)

    def get_channel_id(self, idx):
        try:
            return self.m_vecpszChannelIdList[idx]
        except Exception:
            return None

    def get_channel_name(self, idx):
        return None

    def get_channel_depth(self, idx):
        return None

    def get_channel_bit_count(self, idx):
        return None

    def get_channel_LUT(self, chIdx, lutIdx):
        return None

    def print(self):
        print("Channel Information")
        for cnt in range(len(self.m_vecpszChannelIdList)):
            print(f"\tID = {self.m_vecpszChannelIdList[cnt]}")
            print(f"\tName = {self.m_vecpszChannelName[cnt]}")
            print(f"\tDepth[byte] = {self.m_vecnChannelDepth[cnt]}")
            print(f"\tAvai Bit[bit] = {self.m_vecnChannelBitCount[cnt]}")

    def get_values(self):
        results = []
        for cnt in range(len(self.m_vecpszChannelIdList)):
            result = {
                "id": self.m_vecpszChannelIdList[cnt],
                "name": self.m_vecpszChannelName[cnt],
                "depth": self.m_vecnChannelDepth[cnt],
                "bit_count": self.m_vecnChannelBitCount[cnt],
            }
            results.append(result)
        return results

    @property
    def depth_of_ch0(self):
        return self.m_vecnChannelDepth[0]
