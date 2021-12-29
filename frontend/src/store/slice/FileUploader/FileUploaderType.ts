export const FILE_UPLOADER_SLICE_NAME = 'fileUploader'

export type FileUploader = {
  [id: string]: FileUploaderType
}

export type FileUploaderType = {
  path?: string
  fileName?: string
  isUninitialized: boolean
  pending: boolean
  fulfilled: boolean
  uploadingProgress?: number
  error?: string
}
