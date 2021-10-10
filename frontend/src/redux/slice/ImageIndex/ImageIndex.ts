import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { NODE_DATA_TYPE_SET } from 'const/NodeData'
import { INITIAL_IMAGE_ELEMENT_ID } from 'const/flowchart'
import { uploadImageFile } from './ImageIndexAction'
import { ImageIndex, IMAGE_INDEX_SLICE_NAME } from './ImageIndexType'
import { clickNode } from '../Element/ElementAction'

const initialState: ImageIndex = {
  currentImageId: INITIAL_IMAGE_ELEMENT_ID,
  index: {},
}

export const imageIndexSlice = createSlice({
  name: IMAGE_INDEX_SLICE_NAME,
  initialState,
  reducers: {
    incrementPageIndex: (
      state,
      action: PayloadAction<{ id: string; amount?: number }>,
    ) => {
      const { id, amount } = action.payload
      const newIndex = state.index[id].pageIndex + (amount ?? 1)
      if (newIndex <= state.index[id].maxIndex) {
        state.index[id].pageIndex = newIndex
      }
    },
    decrementPageIndex: (
      state,
      action: PayloadAction<{ id: string; amount?: number }>,
    ) => {
      const { id, amount } = action.payload
      const newIndex = state.index[id].pageIndex - (amount ?? 1)
      if (newIndex >= 0) {
        state.index[id].pageIndex = newIndex
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(uploadImageFile.pending, (state, action) => {
        const { elementId, fileName } = action.meta.arg
        state.currentImageId = elementId
        state.index[elementId] = {
          fileName,
          maxIndex: 0,
          folder: '',
          pageIndex: 0,
          isFulfilled: false,
        }
      })
      .addCase(uploadImageFile.fulfilled, (state, action) => {
        const { elementId, fileName } = action.meta.arg
        const { folderName: folder, maxIndex } = action.payload
        state.currentImageId = elementId
        state.index[elementId] = {
          fileName,
          maxIndex,
          folder,
          pageIndex: 0,
          isFulfilled: true,
        }
      })
      .addCase(clickNode, (state, action) => {
        if (action.payload.type === NODE_DATA_TYPE_SET.DATA) {
          state.currentImageId = action.payload.id
        }
      })
  },
})

export const { incrementPageIndex, decrementPageIndex } =
  imageIndexSlice.actions

export default imageIndexSlice.reducer
