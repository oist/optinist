import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { FileUploader, FILE_UPLOADER_SLICE_NAME } from './FileUploaderType'
import { setUploadProgress, uploadFile } from './FileUploaderActions'
import { inistialUploaderState } from './FileUploaderInitlalState'

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
        const { progess, requestId } = action.payload
        state[requestId].uploadingProgress = progess
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
  },
})

export const { setFileUploaderStateById } = fileUploaderSlice.actions

export default fileUploaderSlice.reducer
