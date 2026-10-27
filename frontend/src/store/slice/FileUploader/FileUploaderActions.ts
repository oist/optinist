import { AxiosProgressEvent } from "axios"

import { createAsyncThunk, createAction } from "@reduxjs/toolkit"

import {
  getStatusLoadViaUrlApi,
  GetStatusViaUrl,
  updateShapeApi,
  uploadFileApi,
  uploadViaUrlApi,
  FILE_TREE_TYPE_SET,
} from "api/files/Files"
import { getFilesTree } from "store/slice/FilesTree/FilesTreeAction"
import { FILE_UPLOADER_SLICE_NAME } from "store/slice/FileUploader/FileUploaderType"
import { FILE_TYPE } from "store/slice/InputNode/InputNodeType"

export const setUploadProgress = createAction<{
  requestId: string
  progress: number
  total: number
}>(`${FILE_UPLOADER_SLICE_NAME}/setUploadProgress`)

export const setDownloadProgress = createAction<{
  requestId: string
  progress: number
}>(`${FILE_UPLOADER_SLICE_NAME}/setDownloadProgress`)

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

export const updateShape = createAsyncThunk<
  boolean,
  {
    workspaceId: number
    fileName: string
  }
>(
  `${FILE_UPLOADER_SLICE_NAME}/updateShape`,
  async ({ workspaceId, fileName }, thunkAPI) => {
    const { dispatch } = thunkAPI
    try {
      const response = await updateShapeApi(workspaceId, fileName)
      dispatch(
        getFilesTree({ workspaceId, fileType: FILE_TREE_TYPE_SET.IMAGE }),
      )
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const uploadViaUrl = createAsyncThunk<
  { file_name: string },
  { workspaceId: number; url: string; requestId: string }
>(`${FILE_UPLOADER_SLICE_NAME}/uploadFileViaUrl`, async (data, thunkAPI) => {
  const { workspaceId, url } = data
  try {
    const data = await uploadViaUrlApi(workspaceId, url)
    return data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getStatusLoadViaUrl = createAsyncThunk<
  GetStatusViaUrl,
  { workspaceId: number; file_name: string; requestId: string }
>(`${FILE_UPLOADER_SLICE_NAME}/getStatusLoadViaUrl`, async (data, thunkAPI) => {
  const { workspaceId, file_name, requestId } = data
  try {
    const data = await getStatusLoadViaUrlApi(workspaceId, file_name)
    thunkAPI.dispatch(
      setUploadProgress({
        requestId,
        progress: (data.current / data.total) * 100,
        total: data.total,
      }),
    )
    return data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
