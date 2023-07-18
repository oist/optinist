import { isNode, Node } from 'react-flow-renderer'
import {
  AlgorithmNodeData,
  NodeData,
  NODE_TYPE_SET,
  InputNodeData,
} from './FlowElementType'

export function isNodeData(
  node: Node<NodeData> | undefined,
): node is Node<NodeData> {
  return node != null && isNode(node) && node.data != null
}

export function isAlgorithmNodeData(
  node: Node<NodeData> | undefined,
): node is Node<AlgorithmNodeData> {
  return isNodeData(node) && node.data?.type === NODE_TYPE_SET.ALGORITHM
}

export function isInputNodeData(
  node: Node<NodeData> | undefined,
): node is Node<InputNodeData> {
  return isNodeData(node) && node.data?.type === NODE_TYPE_SET.INPUT
}

export function getLabelByPath(filePath: string | string[]) {
  if (Array.isArray(filePath)) {
    if (filePath.length === 0) {
      return ''
    } else if (filePath.length === 1) {
      return getFileName(filePath[0])
    } else {
      return getFileName(filePath[0]) + ` ... and ${filePath.length - 1} files`
    }
  } else {
    return getFileName(filePath)
  }
}

export function getFileName(filePath: string) {
  return filePath.split('/').reverse()[0]
}
