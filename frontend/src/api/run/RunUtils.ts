import type { Node } from "reactflow"

import type {
  NodePostDataType,
  InputNodePostData,
  AlgorithmNodePostData,
} from "api/run/Run"
import {
  isInputNodeData,
  isAlgorithmNodeData,
} from "store/slice/FlowElement/FlowElementUtils"

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
