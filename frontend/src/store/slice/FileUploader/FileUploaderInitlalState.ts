import { FileUploaderType } from './FileUploaderType'

export const inistialUploaderState: FileUploaderType = {
  path: undefined,
  fileName: undefined,
  isUninitialized: true,
  pending: false,
  fulfilled: false,
  uploadingProgress: undefined,
  error: undefined,
}
