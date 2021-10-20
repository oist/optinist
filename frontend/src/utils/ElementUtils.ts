import { FlowElement, isNode, Node } from 'react-flow-renderer'
import {
  NodeData,
  InputNodeData,
  AlgoNodeData,
  // OutPutNodeData,
} from 'const/NodeData'

export function isNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<NodeData> {
  return node != null && isNode(node) && node.data != null
}

export function isInputNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<InputNodeData> {
  return isNodeData(node) && node.data != null && node.data.type === 'data'
}

export function isAlgoNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<AlgoNodeData> {
  return isNodeData(node) && node.data != null && node.data.type === 'algo'
}

// export function isOutPutNodeData(
//   node: FlowElement<NodeData> | undefined,
// ): node is Node<OutPutNodeData> {
//   return isNodeData(node) && node.data != null && node.data.type === 'output'
// }
