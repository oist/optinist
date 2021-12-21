export const INPUT_NODE_SLICE_NAME = 'inputNode'

export const FILE_TYPE_SET = {
  CSV: 'csv',
  IMAGE: 'image',
  // NWB:"nwb"
  // JSON:"json"
} as const

export type FILE_TYPE = typeof FILE_TYPE_SET[keyof typeof FILE_TYPE_SET]

export type InputNode = {
  [nodeId: string]: InputNodeType
}

export type InputNodeType = CsvInputNode | ImageInputNode

interface InputNodeBaseType<T extends FILE_TYPE> {
  fileType: T
  selectedFilePath?: string // 複数ファイル指定の予定あり
}

export interface CsvInputNode extends InputNodeBaseType<'csv'> {}

export interface ImageInputNode extends InputNodeBaseType<'image'> {
  maxIndex: number
}
