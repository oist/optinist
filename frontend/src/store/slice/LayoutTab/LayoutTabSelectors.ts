import { RootState } from 'store/store'
import { DisplayTab, ParamFormTab } from './LayoutTabType'

import { isDisplayTab, isParamFormTab } from './LayputTabUtils'

export const selectLayoutTabs = (state: RootState) => state.layoutTab.tabs

export const selectLayoutTab = (tabId: string) => (state: RootState) =>
  selectLayoutTabs(state)[tabId]

export const selectNodeIdByTabId = (tabId: string) => (state: RootState) =>
  selectLayoutTab(tabId)(state).nodeId

export const selectTabTypeByTabId = (tabId: string) => (state: RootState) =>
  selectLayoutTab(tabId)(state).type

export const selectFilePathByTabId = (tabId: string) => (state: RootState) => {
  const tab = selectLayoutTab(tabId)(state)
  if (tab.type === 'displayData') {
    return tab.filePath
  } else {
    return undefined
  }
}
export const selectDataTypeByTabId = (tabId: string) => (state: RootState) => {
  const tab = selectLayoutTab(tabId)(state)
  if (tab.type === 'displayData') {
    return tab.dataType
  } else {
    return undefined
  }
}

export const selectTabByNodeId = (nodeId: string) => (state: RootState) =>
  Object.values(selectLayoutTabs(state)).find((tab) => tab.nodeId === nodeId)

export const selectTabIdByNodeId = (nodeId: string) => (state: RootState) =>
  Object.values(selectLayoutTabs(state)).find((tab) => tab.nodeId === nodeId)
    ?.tabId

export const selectTabIdAndNodeIdList = (state: RootState) =>
  Object.values(selectLayoutTabs(state)).map((tab) => ({
    tabId: tab.tabId,
    nodeId: tab.nodeId,
  }))

export const selectDisplayTabList = (state: RootState) => {
  const displayTabs: DisplayTab[] = []
  for (const tab of Object.values(selectLayoutTabs(state))) {
    if (isDisplayTab(tab)) {
      displayTabs.push(tab)
    }
  }
  return displayTabs
}

export const selectParamFormTabList = (state: RootState) => {
  const paramFormTabs: ParamFormTab[] = []
  for (const tab of Object.values(selectLayoutTabs(state))) {
    if (isParamFormTab(tab)) {
      paramFormTabs.push(tab)
    }
  }
  return paramFormTabs
}

export const selectTabIsList = (state: RootState) =>
  Object.keys(selectLayoutTabs(state))
