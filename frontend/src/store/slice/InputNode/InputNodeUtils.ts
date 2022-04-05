import {
  CsvInputNode,
  ImageInputNode,
  HDF5InputNode,
  InputNodeType,
  FILE_TYPE_SET,
  FluoInputNode,
} from './InputNodeType'

export function isImageInputNode(
  inputNode: InputNodeType,
): inputNode is ImageInputNode {
  return inputNode.fileType === FILE_TYPE_SET.IMAGE
}

export function isCsvInputNode(
  inputNode: InputNodeType,
): inputNode is CsvInputNode {
  return inputNode.fileType === FILE_TYPE_SET.CSV
}

export function isHDF5InputNode(
  inputNode: InputNodeType,
): inputNode is HDF5InputNode {
  return inputNode.fileType === FILE_TYPE_SET.HDF5
}

export function isFluoInputNode(
  inputNode: InputNodeType,
): inputNode is FluoInputNode {
  return inputNode.fileType === FILE_TYPE_SET.FLUO
}
