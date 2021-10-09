import { createAction } from '@reduxjs/toolkit'

import { IMAGE_INDEX_SLICE_NAME } from './ImageIndexType'

export const uploadImageFile = createAction<{
  elementId: string
  fileName: string
  folder: string
  maxIndex: number
}>(`${IMAGE_INDEX_SLICE_NAME}/uploadImageFile`)
