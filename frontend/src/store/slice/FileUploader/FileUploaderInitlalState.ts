import { FileUploaderType } from "store/slice/FileUploader/FileUploaderType"

export const inistialUploaderState: FileUploaderType = {
  path: undefined,
  fileName: undefined,
  isUninitialized: true,
  pending: false,
  fulfilled: false,
  uploadingProgress: undefined,
  error: undefined,
}
