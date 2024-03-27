from studio.app.common.dataclass.image import ImageData
from studio.app.optinist.dataclass.microscope import MicroscopeData


def microscope_to_img(
    microscope: MicroscopeData, output_dir: str, params: dict = None, **kwargs
) -> dict(microscope_image=ImageData):
    reader = microscope.reader
    raw_stack = reader.get_image_stacks()  # (ch, t, y, x) or (ch, t, z, y, x)

    ch = params.get("ch", 0)
    image = ImageData(
        raw_stack[ch],
        output_dir=output_dir,
        file_name="microscope_image",
    )
    microscope.set_data(image.data)

    return {"microscope_image": image}
