import { createAsyncThunk, createAction } from '@reduxjs/toolkit'
import { uploadFileApi } from 'api/files/Files'
import { FILE_TYPE } from '../InputNode/InputNodeType'

import { FILE_UPLOADER_SLICE_NAME } from './FileUploaderType'

export const setUploadProgress = createAction<{
  requestId: string
  progess: number
  total: number
}>(`${FILE_UPLOADER_SLICE_NAME}/setUploadProgress`)

export const uploadFile = createAsyncThunk<
  {
    resultPath: string
  },
  {
    workspaceId: number
    requestId: string
    nodeId?: string
    fileName: string
    formData: FormData
    fileType?: FILE_TYPE
  }
>(
  `${FILE_UPLOADER_SLICE_NAME}/uploadFile`,
  async ({ workspaceId, requestId, fileName, formData }, thunkAPI) => {
    try {
      const config = getUploadConfig((percent, total) => {
        thunkAPI.dispatch(
          setUploadProgress({
            requestId,
            progess: percent,
            total,
          }),
        )
      })
      const response = await uploadFileApi(
        workspaceId,
        fileName,
        config,
        formData,
      )
      return {
        resultPath: response.file_path,
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
