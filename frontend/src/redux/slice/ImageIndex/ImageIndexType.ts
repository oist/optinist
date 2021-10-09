export const IMAGE_INDEX_SLICE_NAME = 'imageIndex'

export interface ImageIndex {
  currentImageId: string
  index: {
    [id: string]: ImageIndexType
  }
}

export type ImageIndexType = {
  fileName: string
  folder: string
  maxIndex: number
  pageIndex: number
}
