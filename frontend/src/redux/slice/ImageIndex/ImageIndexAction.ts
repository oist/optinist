import { createAsyncThunk } from '@reduxjs/toolkit'

import { IMAGE_INDEX_SLICE_NAME } from './ImageIndexType'

export const uploadImageFile = createAsyncThunk<
  { folderName: string; maxIndex: number },
  { elementId: string; fileName: string; formData: FormData }
>(
  `${IMAGE_INDEX_SLICE_NAME}/uploadImageFile_`,
  async ({ elementId, fileName, formData }, thunkAPI) => {
    try {
      const uploadFolderName = `${fileName}(${elementId})`
      const response = await fetch(
        `http://localhost:8000/upload/${uploadFolderName}`,
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
