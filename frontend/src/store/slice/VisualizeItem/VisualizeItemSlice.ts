import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DATA_TYPE, DATA_TYPE_SET } from '../DisplayData/DisplayDataType'

import {
  DefaultSetItem,
  HeatMapItem,
  ImageItem,
  CsvItem,
  TimeSeriesItem,
  RoiItem,
  ScatterItem,
  VisualaizeItem,
  VISUALIZE_ITEM_TYPE_SET,
  BarItem,
  HDF5Item,
} from './VisualizeItemType'
import {
  isDefaultSetItem,
  isDisplayDataItem,
  isHeatMapItem,
  isImageItem,
  isRoiItem,
  isTimeSeriesItem,
  isCsvItem,
  isScatterItem,
} from './VisualizeItemUtils'
import createColormap from 'colormap'

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
  startIndex: 1,
  endIndex: 10,
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
  roiItem: null,
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
  maxIndex: 0,
}
const heatMapItemInitialValue: HeatMapItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.HEAT_MAP,
  showscale: true,
  colors: [
    { rgb: `rgb(0, 0, 255)`, offset: '0' },
    { rgb: `rgb(200, 200, 200)`, offset: '0.5' },
    { rgb: `rgb(255, 0, 0)`, offset: '1.0' },
  ],
}
const csvItemInitialValue: CsvItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.CSV,
  setColumn: null,
  setIndex: false,
  transpose: false,
}
const roiItemInitialValue: RoiItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.ROI,
  colors: createColormap({
    colormap: 'jet',
    nshades: 10,
    format: 'hex',
    alpha: 1,
  }).map((v, idx) => {
    return { rgb: v, offset: String(idx / 9) }
  }),
}
const scatterItemInitialValue: ScatterItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.SCATTER,
  xIndex: 0,
  yIndex: 1,
}
const barItemInitialValue: BarItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.BAR,
}

const hdf5ItemInitialValue: HDF5Item = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.HDF5,
}

function getDisplayDataItemInitialValue(dataType: DATA_TYPE) {
  switch (dataType) {
    case DATA_TYPE_SET.IMAGE:
      return imageItemInitialValue
    case DATA_TYPE_SET.HEAT_MAP:
      return heatMapItemInitialValue
    case DATA_TYPE_SET.TIME_SERIES:
      return timeSeriesItemInitialValue
    case DATA_TYPE_SET.CSV:
      return csvItemInitialValue
    case DATA_TYPE_SET.ROI:
      return roiItemInitialValue
    case DATA_TYPE_SET.SCATTER:
      return scatterItemInitialValue
    case DATA_TYPE_SET.BAR:
      return barItemInitialValue
    case DATA_TYPE_SET.HDF5:
      return hdf5ItemInitialValue
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
    deleteVisualizeItem: (state, action: PayloadAction<number>) => {
      const itemId = action.payload
      delete state.items[itemId]
      if (itemId === state.selectedItemId) {
        state.selectedItemId = null
      }
    },
    addInitialItem: (state) => {
      const nextId = getMaxItemId(state) + 1
      state.items[nextId] = getDisplayDataItemInitialValue(DATA_TYPE_SET.IMAGE)
      state.selectedItemId = nextId
    },
    selectItem: (state, action: PayloadAction<number>) => {
      state.selectedItemId = action.payload
    },
    setRoiItemFilePath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string | null
      }>,
    ) => {
      const { itemId, filePath, nodeId } = action.payload
      const targetItem = state.items[itemId]
      if (isDefaultSetItem(targetItem)) {
        if (targetItem.imageItem.roiItem != null) {
          targetItem.imageItem.roiItem.filePath = filePath
          targetItem.imageItem.roiItem.nodeId = nodeId
        } else {
          targetItem.imageItem.roiItem = {
            ...roiItemInitialValue,
            filePath,
            nodeId,
          }
        }
      } else if (isImageItem(targetItem)) {
        if (targetItem.roiItem != null) {
          targetItem.roiItem.filePath = filePath
          targetItem.roiItem.nodeId = nodeId
        } else {
          targetItem.roiItem = {
            ...roiItemInitialValue,
            filePath,
            nodeId,
          }
        }
      }
    },
    setFilePath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string | null
        dataType?: string
      }>,
    ) => {
      const { itemId, filePath, nodeId, dataType } = action.payload
      const targetItem = state.items[itemId]
      if (isDefaultSetItem(targetItem)) {
        if (dataType === DATA_TYPE_SET.IMAGE) {
          targetItem.imageItem.filePath = filePath
          targetItem.imageItem.nodeId = nodeId
        } else if (dataType === DATA_TYPE_SET.TIME_SERIES) {
          targetItem.timeSeriesItem.filePath = filePath
          targetItem.timeSeriesItem.nodeId = nodeId
        } else if (dataType === DATA_TYPE_SET.HEAT_MAP) {
          targetItem.heatMapItem.filePath = filePath
          targetItem.heatMapItem.nodeId = nodeId
        }
      } else if (isDisplayDataItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      } else {
        throw new Error('error')
      }
    },
    setImageItemFilePath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string | null
      }>,
    ) => {
      const { itemId, filePath, nodeId } = action.payload
      const targetItem = state.items[itemId]
      if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.filePath = filePath
        targetItem.imageItem.nodeId = nodeId
      } else if (isImageItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      }
    },
    setTimeSeriesItemFilePath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string | null
      }>,
    ) => {
      const { itemId, filePath, nodeId } = action.payload
      const targetItem = state.items[itemId]
      if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.filePath = filePath
        targetItem.timeSeriesItem.nodeId = nodeId
      } else if (isTimeSeriesItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      }
    },
    setHeatMapItemFilePath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string | null
      }>,
    ) => {
      const { itemId, filePath, nodeId } = action.payload
      const targetItem = state.items[itemId]
      if (isDefaultSetItem(targetItem)) {
        targetItem.heatMapItem.filePath = filePath
        targetItem.heatMapItem.nodeId = nodeId
      } else if (isHeatMapItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      }
    },
    setDisplayDataPath: (
      state,
      action: PayloadAction<{
        itemId: number
        filePath: string
        nodeId: string | null
        dataType?: DATA_TYPE
      }>,
    ) => {
      const { itemId, filePath, nodeId, dataType } = action.payload
      const targetItem = state.items[itemId]
      if (isDisplayDataItem(targetItem)) {
        if (dataType != null && targetItem.dataType !== dataType) {
          state.items[itemId] = {
            ...getDisplayDataItemInitialValue(dataType),
            filePath,
            nodeId,
          }
        } else {
          targetItem.filePath = filePath
          targetItem.nodeId = nodeId
        }
      } else {
        throw new Error('invalid VisualaizeItemType')
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
      if (type === VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET) {
        state.items[itemId] = defaultSetItemInitialValue
      } else {
        state.items[itemId] = getDisplayDataItemInitialValue(type)
      }
    },
    toggleItemTypeDefaultSet: (state, action: PayloadAction<number>) => {
      const itemId = action.payload
      if (
        state.items[itemId].itemType === VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
      ) {
        state.items[itemId] = {
          ...getDisplayDataItemInitialValue(DATA_TYPE_SET.IMAGE), // FIXME dataTypeの型をNullableに変更して影響箇所も修正する
        }
      } else {
        state.items[itemId] = defaultSetItemInitialValue
      }
    },
    resetImageActiveIndex: (
      state,
      action: PayloadAction<{ itemId: number }>,
    ) => {
      const { itemId } = action.payload
      const targetItem = state.items[itemId]
      if (isImageItem(targetItem)) {
        targetItem.activeIndex = 0
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.activeIndex = 0
      }
    },
    incrementImageActiveIndex: (
      state,
      action: PayloadAction<{ itemId: number }>,
    ) => {
      const { itemId } = action.payload
      const targetItem = state.items[itemId]
      if (isImageItem(targetItem)) {
        targetItem.activeIndex++
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.activeIndex++
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.activeIndex--
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.showticklabels = action.payload.showticklabels
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.zsmooth = action.payload.zsmooth
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.showline = action.payload.showline
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.showgrid = action.payload.showgrid
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.showscale = action.payload.showscale
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.colors = action.payload.colors
      }
    },
    setImageItemStartIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        startIndex: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.startIndex = action.payload.startIndex
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.startIndex = action.payload.startIndex
      }
    },
    setImageItemEndIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        endIndex: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.endIndex = action.payload.endIndex
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.imageItem.endIndex = action.payload.endIndex
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.offset = action.payload.offset
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.span = action.payload.span
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.showgrid = action.payload.showgrid
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.showline = action.payload.showline
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.showticklabels = action.payload.showticklabels
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.zeroline = action.payload.zeroline
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.xrange.left = action.payload.left
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
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.timeSeriesItem.xrange.right = action.payload.right
      }
    },
    setHeatMapItemShowScale: (
      state,
      action: PayloadAction<{
        itemId: number
        showscale: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isHeatMapItem(targetItem)) {
        targetItem.showscale = action.payload.showscale
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.heatMapItem.showscale = action.payload.showscale
      }
    },
    setHeatMapItemColors: (
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
      if (isHeatMapItem(targetItem)) {
        targetItem.colors = action.payload.colors
      } else if (isDefaultSetItem(targetItem)) {
        targetItem.heatMapItem.colors = action.payload.colors
      }
    },
    setRoiItemColors: (
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
      if (isRoiItem(targetItem)) {
        targetItem.colors = action.payload.colors
      }
    },
    setCsvItemTranspose: (
      state,
      action: PayloadAction<{
        itemId: number
        transpose: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isCsvItem(targetItem)) {
        targetItem.transpose = action.payload.transpose
      }
    },
    setCsvItemSetColumn: (
      state,
      action: PayloadAction<{
        itemId: number
        setColumn: number | null
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isCsvItem(targetItem)) {
        targetItem.setColumn = action.payload.setColumn
      }
    },
    setCsvItemSetIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        setIndex: boolean
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isCsvItem(targetItem)) {
        targetItem.setIndex = action.payload.setIndex
      }
    },
    setScatterItemXIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        xIndex: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isScatterItem(targetItem)) {
        targetItem.xIndex = action.payload.xIndex
      }
    },
    setScatterItemYIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        yIndex: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isScatterItem(targetItem)) {
        targetItem.yIndex = action.payload.yIndex
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
  deleteVisualizeItem,
  selectItem,
  setItemType,
  toggleItemTypeDefaultSet,
  setFilePath,
  setHeatMapItemFilePath,
  setImageItemFilePath,
  setTimeSeriesItemFilePath,
  setRoiItemFilePath,
  setDisplayDataPath,
  resetImageActiveIndex,
  incrementImageActiveIndex,
  decrementImageActiveIndex,
  setImageItemShowticklabels,
  setImageItemZsmooth,
  setImageItemShowLine,
  setImageItemShowGrid,
  setImageItemShowScale,
  setImageItemColors,
  setImageItemStartIndex,
  setImageItemEndIndex,
  setTimeSeriesItemOffset,
  setTimeSeriesItemSpan,
  setTimeSeriesItemShowGrid,
  setTimeSeriesItemShowLine,
  setTimeSeriesItemShowTickLabels,
  setTimeSeriesItemZeroLine,
  setTimeSeriesItemXrangeLeft,
  setTimeSeriesItemXrangeRight,
  setHeatMapItemShowScale,
  setHeatMapItemColors,
  setRoiItemColors,
  setCsvItemTranspose,
  setCsvItemSetColumn,
  setCsvItemSetIndex,
  setScatterItemXIndex,
  setScatterItemYIndex,
} = visualaizeItemSlice.actions

export default visualaizeItemSlice.reducer
