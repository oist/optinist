from fastapi import APIRouter, File, Response, UploadFile, Form
import imageio
import pandas as pd
import os

router = APIRouter()

@router.get("/files")
async  def get_files():
    # todo files以下のディレクトリ構造を返す
    return {""}

@router.post("/files/upload/{fileName}/{inputFileNumer}")
async def create_file(response: Response, fileName: str, element_id: str = Form(...), file: UploadFile = File(...), inputFileNumer: int=1):
    root_dir = os.path.join("files", fileName+"("+element_id+")")
    os.makedirs(root_dir, exist_ok=True)

    contents = await file.read()
    tiff_file_path = os.path.join(root_dir, fileName)
    with open(tiff_file_path, "wb") as f:
        f.write(contents)

    tiffs = imageio.volread(tiff_file_path)[:inputFileNumer]

    images = []
    for i, _img in enumerate(tiffs):
        images.append(_img.tolist())

    json_data_path = os.path.join(root_dir, 'image.json')
    pd.DataFrame(images).to_json(json_data_path, indent=4, orient="values")

    return {"json_data_path": json_data_path, "tiff_file_path": tiff_file_path}
