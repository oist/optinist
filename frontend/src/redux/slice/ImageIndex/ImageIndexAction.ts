import { createAsyncThunk } from '@reduxjs/toolkit'

import { IMAGE_INDEX_SLICE_NAME } from './ImageIndexType'

export const uploadImageFile = createAsyncThunk<
  { pngFolder: string; tiffPath: string; maxIndex: number },
  { elementId: string; fileName: string; formData: FormData }
>(
  `${IMAGE_INDEX_SLICE_NAME}/uploadImageFile`,
  async ({ elementId, fileName, formData }, thunkAPI) => {
    try {
      formData.append('element_id', elementId)
      const response = await fetch(
        `http://localhost:8000/api/upload/${fileName}`,
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
