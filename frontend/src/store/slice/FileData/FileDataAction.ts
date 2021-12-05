import { createAsyncThunk } from '@reduxjs/toolkit'
import { BASE_URL } from 'const/API'

import { FILE_DATA_SLICE_NAME } from './FileDataType'

export const uploadImageFile = createAsyncThunk<
  {
    path: string
  },
  {
    nodeId: string
    fileName: string
    formData: FormData
    maxIndex: number
  }
>(
  `${FILE_DATA_SLICE_NAME}/uploadImageFile`,
  async ({ nodeId, fileName, formData }, thunkAPI) => {
    try {
      formData.append('element_id', nodeId)
      const response = await fetch(`${BASE_URL}/files/upload/${fileName}`, {
        method: 'POST',
        mode: 'cors',
        credentials: 'include',
        body: formData,
      })
      const data = await response.json()
      return {
        path: data.file_path,
      }
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const uploadCsvFile = createAsyncThunk<
  {
    path: string
  },
  {
    nodeId: string
    fileName: string
    formData: FormData
  }
>(
  `${FILE_DATA_SLICE_NAME}/uploadCsvFile`,
  async ({ nodeId, fileName, formData }, thunkAPI) => {
    try {
      formData.append('element_id', nodeId)
      const response = await fetch(`${BASE_URL}/files/upload/${fileName}`, {
        method: 'POST',
        mode: 'cors',
        credentials: 'include',
        body: formData,
      })
      const data = await response.json()
      return {
        path: data.file_path,
      }
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
