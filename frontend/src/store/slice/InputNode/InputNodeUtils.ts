import {
  CsvInputNode,
  ImageInputNode,
  InputNodeType,
  FILE_TYPE_SET,
  NWBInputNode,
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

export function isNwbInputNode(
  inputNode: InputNodeType,
): inputNode is NWBInputNode {
  return inputNode.fileType === FILE_TYPE_SET.NWB
}
