import lib


class ObjectiveLensInfo:
    def __init__(self, hAccessor, hArea):
        # variables
        self.m_szName = None
        self.m_dMagnification = None
        self.m_szImmersion = None
        self.m_dNA = None
        self.m_dWD = None
        self.m_dReflectiveIndex = None

        # Get Objective Lens Name
        result, hProp = lib.get_area_property(hAccessor, hArea, "ObjectiveLensInfo")
        result, pName = lib.get_property_value(hAccessor, hProp, "name")
        if len(pName) > 0:
            self.m_szName = pName[0].value.pszString
            del pName

        # Objective Lens Magnification
        result, pMag = lib.get_property_value(hAccessor, hProp, "magnification")
        if len(pMag) > 0:
            self.m_dMagnification = pMag[0].value.dDouble
            del pMag

        # Objective Lens Immersion
        result, pImmersion = lib.get_property_value(hAccessor, hProp, "immersion")
        if len(pImmersion) > 0:
            self.m_szImmersion = pImmersion[0].value.pszString
            del pImmersion

        # Objective Lens NA
        result, pNA = lib.get_property_value(hAccessor, hProp, "na")
        if len(pNA) > 0:
            self.m_dNA = pNA[0].value.dDouble
            del pNA

        # Objective Lens WD
        result, pWD = lib.get_property_value(hAccessor, hProp, "wd")
        if len(pWD) > 0:
            self.m_dWD = pWD[0].value.dDouble
            del pWD

        # Objective Lens reflective index
        result, pRefraction = lib.get_property_value(hAccessor, hProp, "refraction")
        if len(pRefraction) > 0:
            self.m_dReflectiveIndex = pRefraction[0].value.dDouble
            del pRefraction

        if hProp:
            lib.ida.ReleaseProperty(hAccessor, hProp)

    def print(self):
        print("Objective Lens Information")
        print(f"\tName = {self.m_szName}")
        print(f"\tMagnification = {self.m_dMagnification}")
        print(f"\tImmersion = {self.m_szImmersion}")
        print(f"\tNA = {self.m_dNA}")
        print(f"\tWD = {self.m_dWD}")
        print(f"\tReflective Index = {self.m_dReflectiveIndex}")

    def get_name_tm(self, hAccessor, hArea):
        result, hProp = lib.get_area_property(hAccessor, hArea, "ObjectiveLensInfo")
        result, pName = lib.get_property_value(hAccessor, hProp, "name")
        return pName[0].value.pszString

    def get_immersion_tm(self, hAccessor, hArea):
        result, hProp = lib.get_area_property(hAccessor, hArea, "ObjectiveLensInfo")
        result, pImmersion = lib.get_property_value(hAccessor, hProp, "immersion")
        return pImmersion[0].value.pszString

    def get_magnification_tm(self):
        return self.m_dMagnification

    def get_na_tm(self):
        return self.m_dNA

    def get_wd_tm(self):
        return self.m_dWD

    def get_reflective_index_tm(self):
        return self.m_dReflectiveIndex
