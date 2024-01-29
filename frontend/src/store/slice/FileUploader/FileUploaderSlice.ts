import { createSlice, PayloadAction } from "@reduxjs/toolkit"

import {
  getStatusLoadViaUrl,
  setUploadProgress,
  uploadFile,
} from "store/slice/FileUploader/FileUploaderActions"
import { inistialUploaderState } from "store/slice/FileUploader/FileUploaderInitlalState"
import {
  FileUploader,
  FILE_UPLOADER_SLICE_NAME,
} from "store/slice/FileUploader/FileUploaderType"

const initialState: FileUploader = {}

export const fileUploaderSlice = createSlice({
  name: FILE_UPLOADER_SLICE_NAME,
  initialState,
  reducers: {
    setFileUploaderStateById(state, action: PayloadAction<string>) {
      state[action.payload] = inistialUploaderState
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(setUploadProgress, (state, action) => {
        const { progress, requestId } = action.payload
        state[requestId].uploadingProgress = progress
      })
      .addCase(uploadFile.pending, (state, action) => {
        const { fileName, requestId } = action.meta.arg
        const currentState = state[requestId]
        state[requestId] = {
          ...currentState,
          fileName,
          isUninitialized: false,
          pending: true,
          fulfilled: false,
          uploadingProgress: 0,
        }
      })
      .addCase(uploadFile.fulfilled, (state, action) => {
        const { requestId } = action.meta.arg
        const { resultPath } = action.payload
        const currentState = state[requestId]
        state[requestId] = {
          ...currentState,
          path: resultPath,
          pending: false,
          fulfilled: true,
        }
      })
      .addCase(uploadFile.rejected, (state, action) => {
        const { requestId } = action.meta.arg
        const currentState = state[requestId]
        state[requestId] = {
          ...currentState,
          pending: false,
          fulfilled: false,
          error: action.error.message,
        }
      })
      .addCase(getStatusLoadViaUrl.pending, (state, action) => {
        const { requestId } = action.meta.arg
        const currentState = state[requestId]
        state[requestId] = {
          ...currentState,
          isUninitialized: false,
          pending: true,
          fulfilled: false,
        }
      })
      .addCase(getStatusLoadViaUrl.fulfilled, (state, action) => {
        const { current, total } = action.payload
        const { requestId } = action.meta.arg
        const currentState = state[requestId]
        state[requestId] = {
          ...currentState,
          pending: !(current === total),
          fulfilled: true,
        }
      })
      .addCase(getStatusLoadViaUrl.rejected, (state, action) => {
        const { requestId } = action.meta.arg
        const currentState = state[requestId]
        state[requestId] = {
          ...currentState,
          pending: false,
          fulfilled: false,
          error: action.error.message,
        }
      })
  },
})

export const { setFileUploaderStateById } = fileUploaderSlice.actions

export default fileUploaderSlice.reducer
