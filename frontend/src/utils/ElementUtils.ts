import { FlowElement, isNode, Node } from 'react-flow-renderer'
import {
  NodeData,
  ImageNodeData,
  AlgoNodeData,
  NODE_DATA_TYPE_SET,
  CsvNodeData,
} from 'const/NodeData'

export function isNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<NodeData> {
  return node != null && isNode(node) && node.data != null
}

export function isImageNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<ImageNodeData> {
  return (
    isNodeData(node) &&
    node.data != null &&
    node.data.type === NODE_DATA_TYPE_SET.IMAGE
  )
}

export function isCsvNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<CsvNodeData> {
  return (
    isNodeData(node) &&
    node.data != null &&
    node.data.type === NODE_DATA_TYPE_SET.CSV
  )
}

export function isAlgoNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<AlgoNodeData> {
  return (
    isNodeData(node) &&
    node.data != null &&
    node.data.type === NODE_DATA_TYPE_SET.ALGO
  )
}

// export function isOutPutNodeData(
//   node: FlowElement<NodeData> | undefined,
// ): node is Node<OutPutNodeData> {
//   return isNodeData(node) && node.data != null && node.data.type === 'output'
// }
