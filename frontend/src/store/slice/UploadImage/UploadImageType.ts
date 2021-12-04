export const UPLOAD_IMAGE_SLICE_NAME = 'uploadImage'

export interface UploadImage {
  [id: string]: UploadImageType
}

export type UploadImageType = {
  fileName: string
  jsonPath: string
  maxIndex: number
  isFulfilled: boolean
  isUploading: boolean
}
