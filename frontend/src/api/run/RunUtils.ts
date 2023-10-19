import type { Node } from "reactflow"
import {
  isInputNodeData,
  isAlgorithmNodeData,
} from "store/slice/FlowElement/FlowElementUtils"
import type {
  NodePostDataType,
  InputNodePostData,
  AlgorithmNodePostData,
} from "./Run"

export function isInputNodePostData(
  node: Node<NodePostDataType>,
): node is Node<InputNodePostData> {
  return isInputNodeData(node)
}

export function isAlgorithmNodePostData(
  node: Node<NodePostDataType>,
): node is Node<AlgorithmNodePostData> {
  return isAlgorithmNodeData(node)
}
