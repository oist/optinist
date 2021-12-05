import { createAsyncThunk } from '@reduxjs/toolkit'
import { BASE_URL } from 'const/API'

import { UPLOAD_IMAGE_SLICE_NAME } from './ImageFileType'

export const uploadImageFile = createAsyncThunk<
  {
    tiffFilePath: string
  },
  {
    nodeId: string
    fileName: string
    formData: FormData
    inputFileNumber: number
  }
>(
  `${UPLOAD_IMAGE_SLICE_NAME}/uploadImageFile`,
  async ({ nodeId, fileName, formData, inputFileNumber }, thunkAPI) => {
    try {
      formData.append('element_id', nodeId)
      const response = await fetch(
        `${BASE_URL}/files/upload/${fileName}/${inputFileNumber}`,
        {
          method: 'POST',
          mode: 'cors',
          credentials: 'include',
          body: formData,
        },
      )
      const data = await response.json()
      return {
        tiffFilePath: data.tiff_file_path,
      }
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
