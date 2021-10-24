import { NodeData } from 'const/NodeData'
import { FlowElement } from 'react-flow-renderer'

export function toLayoutTabId(
  nodeId: string,
  componentType: ComponentType,
  name: string,
  ...suffix: string[]
): LayoutTabId {
  return `${nodeId}-${componentType}-${name}${
    suffix.length > 0 ? '-' + suffix.reduce((a, b) => `${a}-${b}`) : ''
  }`
}

export function toLayoutTabIdByNode(
  node: FlowElement<NodeData>,
  componentType: ComponentType,
  ...suffix: string[]
) {
  return toLayoutTabId(
    node.id,
    componentType,
    node.data?.label ?? '',
    ...suffix,
  )
}

type LayoutTabId = `${string}-${ComponentType}-${string}`

export type ComponentType = 'paramForm' | 'output' | 'image'
//   | 'flowchart'
//   | 'sidebar'

export function toLayoutTab(
  nodeId: string,
  component: ComponentType,
  name: string,
  ...suffix: string[]
) {
  return {
    id: toLayoutTabId(nodeId, component, name, ...suffix),
    component,
    name: `${component} ${name}`,
    type: 'tab',
  }
}

export function toLayoutTabByNode(
  node: FlowElement<NodeData>,
  componentType: ComponentType,
  ...suffix: string[]
) {
  return toLayoutTab(node.id, componentType, node.data?.label ?? '', ...suffix)
}
