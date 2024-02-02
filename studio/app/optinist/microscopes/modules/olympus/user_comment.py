"""Olympus IDA wrapper module

* Porting of IDA_Sample/UserComment.h,cpp

"""
import studio.app.optinist.microscopes.modules.olympus.lib as lib


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

    def get_values(self):
        return {
            "comment": self.m_szComment,
        }

    @property
    def comment(self):
        return self.m_szComment
