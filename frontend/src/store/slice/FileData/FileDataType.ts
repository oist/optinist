export const FILE_DATA_SLICE_NAME = 'fileData'

export interface FileData {
  image: {
    [nodeId: string]: ImageFileType
  }
  csv: {
    [node: string]: CsvFileType
  }
}

export type ImageFileType = {
  fileName: string
  path: string
  maxIndex: number
  isFulfilled: boolean
  isUploading: boolean
  progess?: number
}

export type CsvFileType = {
  fileName: string
  path: string
  isFulfilled: boolean
  isUploading: boolean
  progess?: number
}
