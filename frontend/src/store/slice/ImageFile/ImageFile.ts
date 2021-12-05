import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { uploadImageFile } from './ImageFileAction'
import { ImageFile, UPLOAD_IMAGE_SLICE_NAME } from './ImageFileType'

const initialState: ImageFile = {}

export const imageFileSlice = createSlice({
  name: UPLOAD_IMAGE_SLICE_NAME,
  initialState,
  reducers: {
    selectImageFile(
      state,
      action: PayloadAction<{ nodeId: string; path: string; maxIndex: number }>,
    ) {
      const { nodeId, path, maxIndex } = action.payload
      const fileName = getFileNameAndParentDirPath(path)
      state[nodeId] = {
        fileName: fileName,
        maxIndex,
        path: path,
        isFulfilled: true,
        isUploading: false,
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(uploadImageFile.pending, (state, action) => {
        const { nodeId, inputFileNumber } = action.meta.arg
        state[nodeId] = {
          fileName: '',
          maxIndex: inputFileNumber,
          path: '',
          isFulfilled: false,
          isUploading: true,
        }
      })
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { nodeId, fileName, inputFileNumber } = action.meta.arg
        const { tiffFilePath } = action.payload
        state[nodeId] = {
          fileName,
          maxIndex: inputFileNumber,
          path: tiffFilePath,
          isFulfilled: true,
          isUploading: false,
        }
      })
  },
})

function getFileNameAndParentDirPath(filePath: string): string {
  let fileName = filePath
  filePath.split('/')
  const lastIndex = filePath.lastIndexOf('/')
  if (lastIndex !== -1) {
    fileName = filePath.substr(lastIndex + 1)
  }
  return fileName
}

export const { selectImageFile } = imageFileSlice.actions

export default imageFileSlice.reducer
