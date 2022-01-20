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
  isTimeSeriesItem,
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
  maxIndex: 10,
  showticklabels: false,
  showline: true,
  zsmooth: 'best',
  showgrid: false,
  showscale: false,
  colors: [
    { rgb: `rgb(0, 0, 0)`, offset: '0' },
    { rgb: `rgb(255, 255, 255)`, offset: '1.0' },
  ],
  activeIndex: 0,
}
const timeSeriesItemInitialValue: TimeSeriesItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.TIME_SERIES,
  offset: false,
  span: 3,
  showgrid: true,
  showline: true,
  showticklabels: true,
  zeroline: false,
  xrange: {
    left: undefined,
    right: undefined,
  },
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
        nodeId: string | null
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
      if (type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET) {
        state.items[itemId] = getDisplayDataItemInitialValue(type)
      } else {
        state.items[itemId] = defaultSetItemInitialValue
      }

      // if (
      //   isDisplayDataItem(targetItem) &&
      //   type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      // ) {
      //   state.items[itemId] = defaultSetItemInitialValue
      //   targetItem.dataType = type
      //   targetItem.filePath = null
      //   targetItem.nodeId = null
      // }
      // } else if (
      //   isDisplayDataItem(targetItem) &&
      //   type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      // ) {
      //   targetItem.dataType = type
      //   targetItem.filePath = null
      //   targetItem.nodeId = null
      // } else if (
      //   isDefaultSetItem(targetItem) &&
      //   type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      // ) {
      //   state.items[itemId] = getDisplayDataItemInitialValue(type)
      // }
      // if (type !== VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET){
      //   state.items[itemId] = getDisplayDataItemInitialValue(type)
      // }
    },
    incrementImageActiveIndex: (
      state,
      action: PayloadAction<{ itemId: number }>,
    ) => {
      const { itemId } = action.payload
      const targetItem = state.items[itemId]
      if (isImageItem(targetItem)) {
        targetItem.activeIndex++
      }
    },
    decrementImageActiveIndex: (
      state,
      action: PayloadAction<{ itemId: number }>,
    ) => {
      const { itemId } = action.payload
      const targetItem = state.items[itemId]
      if (isImageItem(targetItem)) {
        targetItem.activeIndex--
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
    setImageItemZsmooth: (
      state,
      action: PayloadAction<{
        itemId: number
        zsmooth: string | boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.zsmooth = action.payload.zsmooth
      }
    },
    setImageItemShowLine: (
      state,
      action: PayloadAction<{
        itemId: number
        showline: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.showline = action.payload.showline
      }
    },
    setImageItemShowGrid: (
      state,
      action: PayloadAction<{
        itemId: number
        showgrid: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.showgrid = action.payload.showgrid
      }
    },
    setImageItemShowScale: (
      state,
      action: PayloadAction<{
        itemId: number
        showscale: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.showscale = action.payload.showscale
      }
    },
    setImageItemColors: (
      state,
      action: PayloadAction<{
        itemId: number
        colors: {
          rgb: string
          offset: string
        }[]
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.colors = action.payload.colors
      }
    },
    setImageItemMaxIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        maxIndex: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.maxIndex = action.payload.maxIndex
      }
    },
    setTimeSeriesItemOffset: (
      state,
      action: PayloadAction<{
        itemId: number
        offset: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.offset = action.payload.offset
      }
    },
    setTimeSeriesItemSpan: (
      state,
      action: PayloadAction<{
        itemId: number
        span: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.span = action.payload.span
      }
    },
    setTimeSeriesItemShowGrid: (
      state,
      action: PayloadAction<{
        itemId: number
        showgrid: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.showgrid = action.payload.showgrid
      }
    },
    setTimeSeriesItemShowLine: (
      state,
      action: PayloadAction<{
        itemId: number
        showline: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.showline = action.payload.showline
      }
    },
    setTimeSeriesItemShowTickLabels: (
      state,
      action: PayloadAction<{
        itemId: number
        showticklabels: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.showticklabels = action.payload.showticklabels
      }
    },
    setTimeSeriesItemZeroLine: (
      state,
      action: PayloadAction<{
        itemId: number
        zeroline: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.zeroline = action.payload.zeroline
      }
    },
    setTimeSeriesItemXrangeLeft: (
      state,
      action: PayloadAction<{
        itemId: number
        left: number | undefined
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.xrange.left = action.payload.left
      }
    },
    setTimeSeriesItemXrangeRight: (
      state,
      action: PayloadAction<{
        itemId: number
        right: number | undefined
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.xrange.right = action.payload.right
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
  incrementImageActiveIndex,
  decrementImageActiveIndex,
  setImageItemShowticklabels,
  setImageItemZsmooth,
  setImageItemShowLine,
  setImageItemShowGrid,
  setImageItemShowScale,
  setImageItemColors,
  setImageItemMaxIndex,
  setTimeSeriesItemOffset,
  setTimeSeriesItemSpan,
  setTimeSeriesItemShowGrid,
  setTimeSeriesItemShowLine,
  setTimeSeriesItemShowTickLabels,
  setTimeSeriesItemZeroLine,
  setTimeSeriesItemXrangeLeft,
  setTimeSeriesItemXrangeRight,
} = visualaizeItemSlice.actions

export default visualaizeItemSlice.reducer
