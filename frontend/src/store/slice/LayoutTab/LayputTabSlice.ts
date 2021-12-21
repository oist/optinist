import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DATA_TYPE } from '../DisplayData/DisplayDataType'
import { LAYOUT_TAB_SLICE_NAME, LayoutTab } from './LayoutTabType'

const initialState: LayoutTab = { tabs: {} }

export const flexLayoutTabSlice = createSlice({
  name: LAYOUT_TAB_SLICE_NAME,
  initialState,
  reducers: {
    addDisplayTab: (
      state,
      action: PayloadAction<{
        tabId: string // createActonにしても良いかも
        nodeId: string
        filePath: string
        dataType: DATA_TYPE
      }>,
    ) => {
      const { tabId, nodeId, filePath, dataType } = action.payload
      state.tabs[tabId] = {
        tabId,
        type: 'displayData',
        nodeId,
        filePath,
        dataType,
      }
    },
    addParamFormTab: (
      state,
      action: PayloadAction<{
        tabId: string
        nodeId: string
      }>,
    ) => {
      const { tabId, nodeId } = action.payload
      state.tabs[tabId] = {
        tabId,
        type: 'paramForm',
        nodeId,
      }
    },
    deleteAllTabsByIdList: (state, action: PayloadAction<string[]>) => {
      action.payload.forEach((targetTabId) => {
        delete state.tabs[targetTabId]
      })
    },
    deleteAllTabsByNodeId: (state, action: PayloadAction<string>) => {
      const targetNodeId = action.payload
      const targetTabIdList = Object.entries(state.tabs)
        .filter(([tabId, tab]) => tab.nodeId === targetNodeId)
        .map(([tabId, _]) => tabId)
      targetTabIdList.forEach((targetTabId) => {
        delete state.tabs[targetTabId]
      })
    },
    deleteAllTabsByNodeIdList: (state, action: PayloadAction<string[]>) => {
      const targetNodeIdList = action.payload
      const targetTabIdList = Object.entries(state.tabs)
        .filter(([tabId, tab]) => targetNodeIdList.includes(tab.nodeId))
        .map(([tabId, _]) => tabId)
      targetTabIdList.forEach((targetTabId) => {
        delete state.tabs[targetTabId]
      })
    },
  },
})

export const {
  addDisplayTab,
  addParamFormTab,
  deleteAllTabsByIdList,
  deleteAllTabsByNodeId,
  deleteAllTabsByNodeIdList,
} = flexLayoutTabSlice.actions

export default flexLayoutTabSlice.reducer
