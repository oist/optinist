import { NodeData } from 'const/NodeData'
import { FlowElement } from 'react-flow-renderer'

export function toLayoutTabId(
  nodeId: string,
  componentType: ComponentType,
  name: string,
  suffix?: string,
): LayoutTabId {
  return `${nodeId}${SEPARATOR}${componentType}${SEPARATOR}${name}${
    suffix != null ? `${SEPARATOR}${suffix}` : ''
  }`
}

export function toLayoutTabIdByNode(
  node: FlowElement<NodeData>,
  componentType: ComponentType,
  suffix?: string,
) {
  return toLayoutTabId(node.id, componentType, node.data?.label ?? '', suffix)
}

const SEPARATOR = '--'

type LayoutTabId =
  | `${string}${typeof SEPARATOR}${ComponentType}${typeof SEPARATOR}${string}`
  | `${string}${typeof SEPARATOR}${ComponentType}${typeof SEPARATOR}${string}${typeof SEPARATOR}${string}`

export type ComponentType = 'paramForm' | 'output' | 'image'
//   | 'flowchart'
//   | 'sidebar'

export function toLayoutTab(
  nodeId: string,
  component: ComponentType,
  name: string,
  suffix?: string,
) {
  return {
    id: toLayoutTabId(nodeId, component, name, suffix),
    component,
    name: `${component} ${name}`,
    type: 'tab',
  }
}

export function toLayoutTabByNode(
  node: FlowElement<NodeData>,
  componentType: ComponentType,
  suffix?: string,
) {
  return toLayoutTab(node.id, componentType, node.data?.label ?? '', suffix)
}

export function getNodeId(layoutTabId: string): string | null {
  const result = layoutTabId.split(SEPARATOR)[0]
  if (result != null) {
    return result
  } else {
    console.log('failed to get nodeId　from layoutTabId. (' + layoutTabId + ')')
    return null
  }
}

export function getSuffix(layoutTabId: string): string | null {
  const result = layoutTabId.split(SEPARATOR)[3]
  if (result != null) {
    return result
  } else {
    console.log('failed to get suffix from layoutTabId. (' + layoutTabId + ')')
    return null
  }
}
