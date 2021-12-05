import { createSlice } from '@reduxjs/toolkit'

import { uploadImageFile } from '../ImageFile/ImageFileAction'
import { getFiles } from './FilesAction'
import { Files, FILES_SLICE_NAME } from './FilesType'
import { convertToTreeNodeType } from './FilesUtils'

const initialState: Files = {
  isLatest: false,
  isLoading: false,
  tree: [],
}
export const filesSlice = createSlice({
  name: FILES_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getFiles.pending, (state) => {
        state.isLoading = true
      })
      .addCase(getFiles.fulfilled, (state, action) => {
        state.tree = convertToTreeNodeType(action.payload)
        state.isLatest = true
        state.isLoading = false
      })
      .addCase(uploadImageFile.fulfilled, (state) => {
        state.isLatest = false
      })
  },
})

export default filesSlice.reducer
