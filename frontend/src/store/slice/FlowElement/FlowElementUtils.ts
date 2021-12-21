import { FlowElement, isNode, Node } from 'react-flow-renderer'
import { NodeData } from './FlowElementType'

export function isNodeData(
  node: FlowElement<NodeData> | undefined,
): node is Node<NodeData> {
  return node != null && isNode(node) && node.data != null
}
