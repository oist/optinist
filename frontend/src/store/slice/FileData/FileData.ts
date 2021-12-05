import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { uploadCsvFile, uploadImageFile } from './FileDataAction'
import { FileData, FILE_DATA_SLICE_NAME } from './FileDataType'

const initialState: FileData = {
  image: {},
  csv: {},
}

export const fileDataSlice = createSlice({
  name: FILE_DATA_SLICE_NAME,
  initialState,
  reducers: {
    selectImageFile(
      state,
      action: PayloadAction<{ nodeId: string; path: string; maxIndex: number }>,
    ) {
      const { nodeId, path, maxIndex } = action.payload
      const fileName = getFileNameAndParentDirPath(path)
      state.image[nodeId] = {
        fileName: fileName,
        maxIndex,
        path: path,
        isFulfilled: true,
        isUploading: false,
      }
    },
    selectCsvFile(
      state,
      action: PayloadAction<{ nodeId: string; path: string }>,
    ) {
      const { nodeId, path } = action.payload
      const fileName = getFileNameAndParentDirPath(path)
      state.csv[nodeId] = {
        fileName: fileName,
        path: path,
        isFulfilled: true,
        isUploading: false,
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(uploadImageFile.pending, (state, action) => {
        const { nodeId, maxIndex } = action.meta.arg
        state.image[nodeId] = {
          fileName: '',
          maxIndex,
          path: '',
          isFulfilled: false,
          isUploading: true,
        }
      })
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { nodeId, fileName, maxIndex } = action.meta.arg
        const { path } = action.payload
        state.image[nodeId] = {
          fileName,
          maxIndex,
          path,
          isFulfilled: true,
          isUploading: false,
        }
      })
      .addCase(uploadCsvFile.pending, (state, action) => {
        const { nodeId } = action.meta.arg
        state.csv[nodeId] = {
          fileName: '',
          path: '',
          isFulfilled: false,
          isUploading: true,
        }
      })
      .addCase(uploadCsvFile.fulfilled, (state, action) => {
        const { nodeId, fileName } = action.meta.arg
        const { path } = action.payload
        state.csv[nodeId] = {
          fileName,
          path,
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

export const { selectImageFile, selectCsvFile } = fileDataSlice.actions

export default fileDataSlice.reducer
