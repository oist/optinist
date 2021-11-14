import { createSlice } from '@reduxjs/toolkit'

import { uploadImageFile } from './UploadImageAction'
import { UploadImage, UPLOAD_IMAGE_SLICE_NAME } from './UploadImageType'

const initialState: UploadImage = {}

export const uploadImageSlice = createSlice({
  name: UPLOAD_IMAGE_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(uploadImageFile.pending, (state, action) => {
        const { nodeId, fileName, inputFileNumber } = action.meta.arg
        state[nodeId] = {
          fileName,
          maxIndex: inputFileNumber,
          path: '',
          isFulfilled: false,
        }
      })
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { nodeId, fileName, inputFileNumber } = action.meta.arg
        const { jsonDataPath } = action.payload
        state[nodeId] = {
          fileName,
          maxIndex: inputFileNumber,
          path: jsonDataPath,
          isFulfilled: true,
        }
      })
  },
})

export default uploadImageSlice.reducer
