import { FlowElement, isNode, Node } from 'react-flow-renderer'
import {
  AlgorithmNodeData,
  NodeData,
  NODE_TYPE_SET,
  InputNodeData,
} from './FlowElementType'

export function isNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<NodeData> {
  return node != null && isNode(node) && node.data != null
}

export function isAlgorithmNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<AlgorithmNodeData> {
  return isNodeData(node) && node.data?.type === NODE_TYPE_SET.ALGORITHM
}

export function isInputNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<InputNodeData> {
  return isNodeData(node) && node.data?.type === NODE_TYPE_SET.INPUT
}
