import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DATA_TYPE, DATA_TYPE_SET } from '../DisplayData/DisplayDataType'

import {
  DefaultSetItem,
  HeatMapItem,
  ImageItem,
  TableItem,
  TimeSeriesItem,
  VisualaizeItem,
  VISUALIZE_ITEM_TYPE_SET,
} from './VisualizeItemType'
import {
  isDefaultSetItem,
  isDisplayDataItem,
  isImageItem,
} from './VisualizeItemUtils'
export const initialState: VisualaizeItem = {
  items: {},
  selectedItemId: null,
}
const displayDataCommonInitialValue = {
  itemType: VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA,
  filePath: null,
  nodeId: null,
}
const imageItemInitialValue: ImageItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.IMAGE,
  showticklabels: false,
}
const timeSeriesItemInitialValue: TimeSeriesItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.TIME_SERIES,
}
const heatMapItemInitialValue: HeatMapItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.HEAT_MAP,
}
const tableItemInitialValue: TableItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.TABLE,
}
function getDisplayDataItemInitialValue(dataType: DATA_TYPE) {
  switch (dataType) {
    case DATA_TYPE_SET.IMAGE:
      return imageItemInitialValue
    case DATA_TYPE_SET.HEAT_MAP:
      return heatMapItemInitialValue
    case DATA_TYPE_SET.TIME_SERIES:
      return timeSeriesItemInitialValue
    case DATA_TYPE_SET.TABLE:
      return tableItemInitialValue
  }
}

const defaultSetItemInitialValue: DefaultSetItem = {
  itemType: VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET,
  imageItem: imageItemInitialValue,
  timeSeriesItem: timeSeriesItemInitialValue,
  heatMapItem: heatMapItemInitialValue,
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
        state.items[itemId] = getDisplayDataItemInitialValue(type)
      }
    },
    setImageItemShowticklabels: (
      state,
      action: PayloadAction<{
        itemId: number
        showticklabels: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.showticklabels = action.payload.showticklabels
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
  setImageItemShowticklabels,
} = visualaizeItemSlice.actions

export default visualaizeItemSlice.reducer