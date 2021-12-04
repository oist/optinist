import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { uploadImageFile } from './UploadImageAction'
import { UploadImage, UPLOAD_IMAGE_SLICE_NAME } from './UploadImageType'

const initialState: UploadImage = {}

export const uploadImageSlice = createSlice({
  name: UPLOAD_IMAGE_SLICE_NAME,
  initialState,
  reducers: {
    selectImageFile(
      state,
      action: PayloadAction<{ nodeId: string; path: string; maxIndex: number }>,
    ) {
      const { nodeId, path, maxIndex } = action.payload
      const [parentDirPath, fileName] = getFileNameAndParentDirPath(path)
      state[nodeId] = {
        fileName: fileName,
        maxIndex,
        jsonPath: parentDirPath + '/image.json',
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
          jsonPath: '',
          isFulfilled: false,
          isUploading: true,
        }
      })
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { nodeId, fileName, inputFileNumber } = action.meta.arg
        const { jsonDataPath } = action.payload
        state[nodeId] = {
          fileName,
          maxIndex: inputFileNumber,
          jsonPath: jsonDataPath,
          isFulfilled: true,
          isUploading: false,
        }
      })
  },
})

function getFileNameAndParentDirPath(filePath: string): [string, string] {
  let dirPath = ''
  let fileName = filePath
  filePath.split('/')
  const lastIndex = filePath.lastIndexOf('/')
  if (lastIndex !== -1) {
    dirPath = filePath.substr(0, lastIndex)
    fileName = filePath.substr(lastIndex + 1)
  }
  return [dirPath, fileName]
}

export const { selectImageFile } = uploadImageSlice.actions

export default uploadImageSlice.reducer
