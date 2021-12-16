import { createAsyncThunk, createAction } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'

import { FileData, FILE_DATA_SLICE_NAME } from './FileDataType'

export const setUploadProgress = createAction<{
  nodeId: string
  fileDataKey: keyof FileData
  progess: number
  total: number
}>(FILE_DATA_SLICE_NAME + '/setUploadProgress ')

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
      const config = getUploadConfig((percent, total) => {
        thunkAPI.dispatch(
          setUploadProgress({
            nodeId,
            fileDataKey: 'image',
            progess: percent,
            total,
          }),
        )
      })
      const response = await axios.post(
        `${BASE_URL}/files/upload/${fileName}`,
        formData,
        config,
      )
      const data = response.data
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
      const config = getUploadConfig((percent, total) => {
        thunkAPI.dispatch(
          setUploadProgress({
            nodeId,
            fileDataKey: 'csv',
            progess: percent,
            total,
          }),
        )
      })
      const response = await axios.post(
        `${BASE_URL}/files/upload/${fileName}`,
        formData,
        config,
      )
      const data = response.data
      return {
        path: data.file_path,
      }
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

function getUploadConfig(
  onUpdateProgressFn: (percent: number, totalSize: number) => void,
) {
  return {
    onUploadProgress: function (progressEvent: any) {
      const percentCompleted = Math.round(
        (progressEvent.loaded * 100) / progressEvent.total,
      )
      onUpdateProgressFn(percentCompleted, progressEvent.total)
    },
  }
}
