import { AxiosProgressEvent } from "axios"

import { createAsyncThunk, createAction } from "@reduxjs/toolkit"

import { uploadFileApi, uploadViaUrlApi } from "api/files/Files"
import { FILE_UPLOADER_SLICE_NAME } from "store/slice/FileUploader/FileUploaderType"
import { FILE_TYPE } from "store/slice/InputNode/InputNodeType"

export const setUploadProgress = createAction<{
  requestId: string
  progress: number
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
            progress: percent,
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
    onUploadProgress: function (progressEvent: AxiosProgressEvent) {
      if (!progressEvent.total) {
        onUpdateProgressFn(0, 100)
      } else {
        const percentCompleted = Math.round(
          (progressEvent.loaded * 100) / progressEvent.total,
        )
        onUpdateProgressFn(percentCompleted, progressEvent.total)
      }
    },
  }
}

export const uploadViaUrl = createAsyncThunk<
  { file_name: string },
  { workspaceId: number; url: string; requestId: string }
>(`${FILE_UPLOADER_SLICE_NAME}/uploadFileViaUrl`, async (data, thunkAPI) => {
  const { workspaceId, url, requestId } = data
  try {
    const config = getUploadConfig((percent, total) => {
      thunkAPI.dispatch(
        setUploadProgress({
          requestId,
          progress: percent,
          total,
        }),
      )
    })
    const data = await uploadViaUrlApi(workspaceId, url, config)
    return data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
