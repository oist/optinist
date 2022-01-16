import { createSlice, PayloadAction, createAction } from '@reduxjs/toolkit'
import { DATA_TYPE, DATA_TYPE_SET } from '../DisplayData/DisplayDataType'

import {
  DefaultSetItem,
  DisplayDataItem,
  VisualaizeItem,
  VISUALIZE_ITEM_TYPE,
  VISUALIZE_ITEM_TYPE_SET,
} from './VisualizeItemType'
import { isDefaultSetItem, isDisplayDataItem } from './VisualizeItemUtils'
export const initialState: VisualaizeItem = {
  items: {},
  selectedItemId: null,
}
const displayDataItemInitialValue: DisplayDataItem = {
  itemType: VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA,
  dataType: null,
  filePath: null,
  nodeId: null,
}
const defaultSetItemInitialValue: DefaultSetItem = {
  itemType: VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET,
  imageItem: {
    itemType: VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA,
    dataType: DATA_TYPE_SET.IMAGE,
    filePath: null,
    nodeId: null,
  },
  timeSeriesItem: {
    itemType: VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA,
    dataType: DATA_TYPE_SET.TIME_SERIES,
    filePath: null,
    nodeId: null,
  },
  heatMapItem: {
    itemType: VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA,
    dataType: DATA_TYPE_SET.HEAT_MAP,
    filePath: null,
    nodeId: null,
  },
}
export const visualaizeItemSlice = createSlice({
  name: 'visualaizeItem',
  initialState,
  reducers: {
    deleteItem: (state, action: PayloadAction<number>) => {
      const itemId = action.payload
      delete state.items[itemId]
      if (itemId === state.selectedItemId) {
        state.selectedItemId = null
      }
    },
    addInitialItem: (state) => {
      const nextId = getMaxItemId(state) + 1
      state.items[nextId] = defaultSetItemInitialValue
      state.selectedItemId = nextId
    },
    selectItem: (state, action: PayloadAction<number>) => {
      state.selectedItemId = action.payload
    },
    setDisplayDataPath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string
      }>,
    ) => {
      const { itemId, filePath, nodeId } = action.payload
      const targetItem = state.items[itemId]
      if (isDisplayDataItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      }
    },
    setItemType: (
      state,
      action: PayloadAction<{
        itemId: number
        type: typeof VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET | DATA_TYPE
      }>,
    ) => {
      const { itemId, type } = action.payload
      const targetItem = state.items[itemId]
      if (
        isDisplayDataItem(targetItem) &&
        type === VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      ) {
        state.items[itemId] = defaultSetItemInitialValue
      } else if (
        isDisplayDataItem(targetItem) &&
        type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      ) {
        targetItem.dataType = type
        targetItem.filePath = null
        targetItem.nodeId = null
      } else if (
        isDefaultSetItem(targetItem) &&
        type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      ) {
        state.items[itemId] = {
          ...displayDataItemInitialValue,
          dataType: type,
        }
      }
    },
  },
})

function getMaxItemId(state: VisualaizeItem) {
  const idList = Object.keys(state.items).map((key) => Number(key))
  const maxId = idList.length > 0 ? idList.reduce((a, b) => Math.max(a, b)) : 0
  return maxId
}

export const {
  addInitialItem,
  deleteItem,
  selectItem,
  setItemType,
  setDisplayDataPath,
} = visualaizeItemSlice.actions

export default visualaizeItemSlice.reducer
