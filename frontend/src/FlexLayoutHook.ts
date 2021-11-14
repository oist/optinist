import React from 'react'
import { useSelector } from 'react-redux'

import { FlexLayoutModelContext } from './App'
import {
  ComponentType,
  toLayoutTabByNode,
  toLayoutTabIdByNode,
} from 'utils/FlexLayoutUtils'
import { Actions, DockLocation } from 'flexlayout-react'
import { nodeByIdSelector } from 'store/slice/Element/ElementSelector'

export function useTabAction(nodeId: string) {
  const node = useSelector(nodeByIdSelector(nodeId))
  const model = React.useContext(FlexLayoutModelContext)
  if (node === undefined) {
    return null
  }
  return (component: ComponentType, toTabSetId: string, suffix?: string) => {
    const layoutId = toLayoutTabIdByNode(node, component, suffix)
    if (model.getNodeById(layoutId) != null) {
      return Actions.selectTab(toLayoutTabIdByNode(node, component, suffix))
    } else {
      return Actions.addNode(
        toLayoutTabByNode(node, component, suffix),
        toTabSetId,
        DockLocation.CENTER,
        0,
      )
    }
  }
}

export function useGetDeleteTabActions() {
  const model = React.useContext(FlexLayoutModelContext)
  return (...nodeIds: string[]) => {
    // layoutの構造に山を貼ってidを取得する findやfilter関数が提供されていればそれが使えるのに...
    const paramIds = model
      .getRoot()
      .getChildren()[0]
      .getChildren()[1]
      .getChildren()[0]
      .getChildren()
      .filter(
        (elm) => nodeIds.filter((id) => elm.getId().startsWith(id)).length > 0,
      )
      .map((elm) => elm.getId())
    const outputIds = model
      .getRoot()
      .getChildren()[0]
      .getChildren()[1]
      .getChildren()[1]
      .getChildren()
      .filter(
        (elm) => nodeIds.filter((id) => elm.getId().startsWith(id)).length > 0,
      )
      .map((elm) => elm.getId())
    const actions = paramIds
      .concat(outputIds)
      .map((id) => Actions.deleteTab(id))
    return actions
  }
}
