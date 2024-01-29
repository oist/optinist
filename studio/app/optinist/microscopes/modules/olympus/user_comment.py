import lib


class UserComment:
    def __init__(self, hAccessor, hArea):
        # variables

        result, hProp = lib.get_area_property(hAccessor, hArea, "UserComment")

        result, pUserComment = lib.get_property_value(hAccessor, hProp, "comment")
        self.m_szComment = pUserComment[0].value.pszString
        del pUserComment

        if hProp:
            lib.ida.ReleaseProperty(hAccessor, hProp)

    def print(self):
        print("User Comment")
        print(f"\tComment = {self.m_szComment}")
