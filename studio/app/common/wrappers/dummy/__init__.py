from studio.app.common.wrappers.dummy.dummy import (
    dummy_image2heat,
    dummy_image2image,
    dummy_image2image8roi,
    dummy_image2image8roi8time8heat,
    dummy_image2image8time,
    dummy_image2roi,
    dummy_image2scatter,
    dummy_image2time,
    dummy_image2time8iscell,
    dummy_keyerror,
    dummy_time2time,
    dummy_typeerror,
)

dummy_wrapper_dict = {
    "dummy": {
        "dummy_image2image": {
            "function": dummy_image2image,
        },
        "dummy_image2time": {
            "function": dummy_image2time,
        },
        "dummy_image2heat": {
            "function": dummy_image2heat,
        },
        "dummy_time2time": {
            "function": dummy_time2time,
        },
        "dummy_image2image8time": {
            "function": dummy_image2image8time,
        },
        "dummy_keyerror": {
            "function": dummy_keyerror,
        },
        "dummy_typeerror": {
            "function": dummy_typeerror,
        },
        "dummy_image2time8iscell": {
            "function": dummy_image2time8iscell,
        },
        "dummy_image2roi": {
            "function": dummy_image2roi,
        },
        "dummy_image2image8roi": {
            "function": dummy_image2image8roi,
        },
        "dummy_image2image8roi8time8heat": {
            "function": dummy_image2image8roi8time8heat,
        },
        "dummy_image2scatter": {
            "function": dummy_image2scatter,
        },
    }
}
