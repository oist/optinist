import React from 'react'
import { nanoid } from '@reduxjs/toolkit'
import { useSelector, useDispatch } from 'react-redux'
import { Actions, DockLocation } from 'flexlayout-react'

import {
  selectDisplayTabList,
  selectParamFormTabList,
  selectTabIdAndNodeIdList,
} from 'store/slice/LayoutTab/LayoutTabSelectors'
import {
  DisplayTab,
  ParamFormTab,
  TAB_COMPONENT_TYPE_SET,
} from 'store/slice/LayoutTab/LayoutTabType'
import { selectAlgorithmName } from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import {
  addDisplayTab,
  addParamFormTab,
  deleteAllTabsByIdList,
} from 'store/slice/LayoutTab/LayputTabSlice'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import { equalsDisplayTabInfo } from 'store/slice/LayoutTab/LayputTabUtils'

import { FlexLayoutModelContext } from 'App'
import { PARAM_FORM_TABSET_ID, DISPLAY_DATA_TABSET_ID } from 'const/flexlayout'
import { toLayoutTab } from './FlexLayputUtils'

export function useDisplayDataTabAciton(nodeId: string) {
  const model = React.useContext(FlexLayoutModelContext)
  const displayTabList = useSelector(
    selectDisplayTabList,
    displayTabListEquality,
  )
  const dispatch = useDispatch()
  return {
    setDisplayTab: React.useCallback(
      (filePath: string, dataType: DATA_TYPE, name?: string) => {
        const sameInfoTab = displayTabList.find((tab) =>
          equalsDisplayTabInfo(tab, { filePath, dataType, nodeId }),
        )
        if (sameInfoTab !== undefined) {
          model.doAction(Actions.selectTab(sameInfoTab.tabId))
        } else {
          const newTabId = nanoid()
          model.doAction(
            Actions.addNode(
              toLayoutTab(
                newTabId,
                TAB_COMPONENT_TYPE_SET.DISPLAY_DATA,
                name ?? getFileName(filePath),
              ),
              DISPLAY_DATA_TABSET_ID,
              DockLocation.CENTER,
              0,
            ),
          )
          dispatch(
            addDisplayTab({ tabId: newTabId, nodeId, filePath, dataType }),
          )
        }
      },
      [model, dispatch, nodeId, displayTabList],
    ),
    deleteDisplayTab: React.useCallback(
      (filePath: string) => {
        const targetTabIdList = displayTabList
          .filter(
            ({ nodeId: tabNodeId, filePath: tabFilePath }) =>
              tabNodeId === nodeId && tabFilePath === filePath,
          )
          .map(({ tabId }) => tabId)
        targetTabIdList.forEach((targetTabId) =>
          model.doAction(Actions.deleteTab(targetTabId)),
        )
        dispatch(deleteAllTabsByIdList(targetTabIdList))
      },
      [model, dispatch, displayTabList, nodeId],
    ),
  }
}

function getFileName(filePath: string) {
  let name = filePath
  const array = filePath.split('/')
  if (array.length > 0) {
    name = array[array.length - 1]
  }
  return name
}

export function useParamFormTabAction(nodeId: string) {
  const model = React.useContext(FlexLayoutModelContext)
  const dispatch = useDispatch()
  const paramFormTabList = useSelector(
    selectParamFormTabList,
    paramFormTabListEquality,
  )
  const name = useSelector(selectAlgorithmName(nodeId))
  return React.useCallback(() => {
    const sameInfoTab = paramFormTabList.find((tab) => tab.nodeId === nodeId)
    if (sameInfoTab !== undefined) {
      model.doAction(Actions.selectTab(sameInfoTab.tabId))
    } else {
      const tabId = nanoid()
      model.doAction(
        Actions.addNode(
          toLayoutTab(tabId, TAB_COMPONENT_TYPE_SET.PARAM_FORM, name),
          PARAM_FORM_TABSET_ID,
          DockLocation.CENTER,
          0,
        ),
      )
      dispatch(addParamFormTab({ tabId, nodeId }))
    }
  }, [model, dispatch, paramFormTabList, nodeId, name])
}

export function useDeleteTabByNodeIdListActions() {
  const model = React.useContext(FlexLayoutModelContext)
  const dispatch = useDispatch()
  const tabIdAndNodeIdList = useSelector(
    selectTabIdAndNodeIdList,
    tabAndNodeIdEquality,
  )
  // targetNodeIdsに紐付いているtabIdがkeyとなっているものを削除
  return (targetNodeIds: string[]) => {
    const targetTabIdList = tabIdAndNodeIdList
      .filter(({ tabId, nodeId }) => targetNodeIds.includes(nodeId))
      .map(({ tabId, nodeId }) => tabId)
    targetTabIdList.forEach((targetTabId) =>
      model.doAction(Actions.deleteTab(targetTabId)),
    )
    dispatch(deleteAllTabsByIdList(targetTabIdList))
  }
}

function displayTabListEquality(a: DisplayTab[], b: DisplayTab[]) {
  return (
    a === b ||
    (a.length === b.length && a.every((v, i) => displayTabEquality(v, b[i])))
  )
}

function displayTabEquality(a: DisplayTab, b: DisplayTab) {
  return (
    a.tabId === b.tabId &&
    a.nodeId === b.nodeId &&
    a.filePath === b.filePath &&
    a.dataType === b.dataType
  )
}

function paramFormTabListEquality(a: ParamFormTab[], b: ParamFormTab[]) {
  return (
    a === b ||
    (a.length === b.length && a.every((v, i) => paramFormTabEquality(v, b[i])))
  )
}

function paramFormTabEquality(a: ParamFormTab, b: ParamFormTab) {
  return a.tabId === b.tabId && a.nodeId === b.nodeId
}

function tabAndNodeIdEquality(
  a: {
    tabId: string
    nodeId: string
  }[],
  b: {
    tabId: string
    nodeId: string
  }[],
) {
  return (
    a === b ||
    (a.length === b.length &&
      a.every((v, i) => v.nodeId === b[i].nodeId && v.tabId === b[i].tabId))
  )
}
