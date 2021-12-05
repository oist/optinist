export const UPLOAD_IMAGE_SLICE_NAME = 'imageFile'

export interface ImageFile {
  [id: string]: ImageFileType
}

export type ImageFileType = {
  fileName: string
  path: string
  maxIndex: number
  isFulfilled: boolean
  isUploading: boolean
}
