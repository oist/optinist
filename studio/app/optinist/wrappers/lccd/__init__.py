from studio.app.optinist.wrappers.lccd.lccd_detection import lccd_detect

lccd_wrapper_dict = {
    "lccd": {
        "lccd_cell_detection": {
            "function": lccd_detect,
            "conda_name": "lccd",
        },
    }
}
