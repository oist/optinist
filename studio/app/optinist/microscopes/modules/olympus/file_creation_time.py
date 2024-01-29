import lib


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

    def get_file_creation_time_tm(self, hAccessor, hArea):
        result, hProp = lib.get_area_property(hAccessor, hArea, "CreationDateTime")
        result, pCreationDateTime = lib.get_property_value(hAccessor, hProp, "dateTime")
        return pCreationDateTime[0].value.pszString
