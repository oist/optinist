import { FlowElement, isNode, Node } from 'react-flow-renderer'
import {
  NodeDataType,
  InputNodeData,
  AlgoNodeData,
  OutPutNodeData,
} from './ElementType'

export function isNodeData(
  node: FlowElement<NodeDataType> | undefined,
): node is Node<NodeDataType> {
  return node != null && isNode(node) && node.data != null
}

export function isInputNodeData(
  node: FlowElement<NodeDataType> | undefined,
): node is Node<InputNodeData> {
  return isNodeData(node) && node.data != null && node.data.type === 'input'
}

export function isAlgoNodeData(
  node: FlowElement<NodeDataType> | undefined,
): node is Node<AlgoNodeData> {
  return isNodeData(node) && node.data != null && node.data.type === 'algo'
}

export function isOutPutNodeData(
  node: FlowElement<NodeDataType> | undefined,
): node is Node<OutPutNodeData> {
  return isNodeData(node) && node.data != null && node.data.type === 'output'
}
