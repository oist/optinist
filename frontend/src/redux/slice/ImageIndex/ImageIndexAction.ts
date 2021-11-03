import { createAsyncThunk } from '@reduxjs/toolkit'
import { BASE_URL } from 'const/API'

import { IMAGE_INDEX_SLICE_NAME } from './ImageIndexType'

export const uploadImageFile = createAsyncThunk<
  { pngFolder: string; tiffPath: string; maxIndex: number },
  {
    elementId: string
    fileName: string
    formData: FormData
    inputFileNumber: number
  }
>(
  `${IMAGE_INDEX_SLICE_NAME}/uploadImageFile`,
  async ({ elementId, fileName, formData, inputFileNumber }, thunkAPI) => {
    try {
      formData.append('element_id', elementId)
      const response = await fetch(
        `${BASE_URL}/api/upload/${fileName}/${inputFileNumber}`,
        {
          method: 'POST',
          mode: 'cors',
          credentials: 'include',
          body: formData,
        },
      )
      const data = await response.json()
      return data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
