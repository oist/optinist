import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DATA_TYPE, DATA_TYPE_SET } from '../DisplayData/DisplayDataType'
import {
  deleteDisplayItem,
  selectingImageArea,
  setImageItemClikedDataId,
  setNewDisplayDataPath,
} from './VisualizeItemActions'

import {
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
  HTMLItem,
  FluoItem,
  BehaviorItem,
  VISUALIZE_ITEM_SLICE_NAME,
} from './VisualizeItemType'
import {
  isDisplayDataItem,
  isHeatMapItem,
  isImageItem,
  isTimeSeriesItem,
  isCsvItem,
  isScatterItem,
  isBarItem,
} from './VisualizeItemUtils'

export const initialState: VisualaizeItem = {
  items: {},
  selectedItemId: null,
  layout: [],
}
const displayDataCommonInitialValue = {
  itemType: VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA,
  filePath: null,
  nodeId: null,
  width: 500,
  height: 500,
  isWorkflowDialog: false,
  saveFileName: 'newPlot',
  saveFormat: 'png',
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
  alpha: 1.0,
  roiItem: null,
  roiAlpha: 1.0,
  duration: 500,
}
const timeSeriesItemInitialValue: TimeSeriesItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.TIME_SERIES,
  offset: true,
  span: 5,
  showgrid: true,
  showline: true,
  showticklabels: true,
  zeroline: false,
  xrange: {
    left: undefined,
    right: undefined,
  },
  maxIndex: 0,
  drawOrderList: [],
  refImageItemId: null,
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
  setHeader: null,
  setIndex: false,
  transpose: false,
}
const roiItemInitialValue: RoiItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.ROI,
}
const scatterItemInitialValue: ScatterItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.SCATTER,
  xIndex: '0',
  yIndex: '1',
}
const barItemInitialValue: BarItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.BAR,
  index: '0',
}
const hdf5ItemInitialValue: HDF5Item = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.HDF5,
}
const htmlItemInitialValue: HTMLItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.HTML,
}
const fluoItemInitialValue: FluoItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.FLUO,
}
const behaviorItemInitialValue: BehaviorItem = {
  ...displayDataCommonInitialValue,
  dataType: DATA_TYPE_SET.BEHAVIOR,
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
    case DATA_TYPE_SET.HTML:
      return htmlItemInitialValue
    case DATA_TYPE_SET.FLUO:
      return fluoItemInitialValue
    case DATA_TYPE_SET.BEHAVIOR:
      return behaviorItemInitialValue
  }
}

export const visualaizeItemSlice = createSlice({
  name: VISUALIZE_ITEM_SLICE_NAME,
  initialState,
  reducers: {
    pushInitialItemToNewRow: (state) => {
      const newItemId = addInitialItemFn(state)
      state.layout.push([newItemId])
    },
    insertInitialItemToNextColumn: (state, action: PayloadAction<number>) => {
      const newItemId = addInitialItemFn(state)
      const targetItemId = action.payload
      const targetRowIndex = state.layout.findIndex((row) =>
        row.includes(targetItemId),
      )
      const targetColumnIndex =
        state.layout[targetRowIndex].indexOf(targetItemId)
      state.layout[targetRowIndex].splice(targetColumnIndex + 1, 0, newItemId)
    },
    addItemForWorkflowDialog: (
      state,
      action: PayloadAction<{
        nodeId: string
        filePath: string
        dataType: DATA_TYPE
      }>,
    ) => {
      const { nodeId, filePath, dataType } = action.payload
      const maxId = getMaxItemId(state)
      const newItemId = maxId != null ? maxId + 1 : 0
      state.items[newItemId] = {
        ...getDisplayDataItemInitialValue(dataType),
        isWorkflowDialog: true,
        nodeId,
        filePath,
      }
    },
    deleteAllItemForWorkflowDialog: (state) => {
      const targetItemIdList = Object.entries(state.items)
        .filter(([itemId, value]) => value.isWorkflowDialog)
        .map(([itemId, value]) => Number(itemId))
      targetItemIdList.forEach(
        (targetItemId) => delete state.items[targetItemId],
      )
    },
    setItemSize: (
      state,
      action: PayloadAction<{
        itemId: number
        width: number
        height: number
      }>,
    ) => {
      const { itemId, width, height } = action.payload
      const targetItem = state.items[itemId]
      targetItem.width = width
      targetItem.height = height
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
        outputKey?: string
      }>,
    ) => {
      const { itemId, filePath, nodeId, outputKey } = action.payload
      const targetItem = state.items[itemId]
      if (isImageItem(targetItem)) {
        Object.values(state.items).forEach((item) => {
          if (
            isTimeSeriesItem(item) &&
            item.filePath != null &&
            item.refImageItemId === itemId
          ) {
            item.drawOrderList = []
          }
        })

        if (targetItem.roiItem != null) {
          targetItem.roiItem.filePath = filePath
          targetItem.roiItem.nodeId = nodeId
          targetItem.roiItem.outputKey = outputKey
        } else {
          targetItem.roiItem = {
            ...roiItemInitialValue,
            filePath,
            nodeId,
            outputKey,
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
      }>,
    ) => {
      const { itemId, filePath, nodeId } = action.payload
      const targetItem = state.items[itemId]
      if (isDisplayDataItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      } else {
        throw new Error('error')
      }
    },
    setSaveFormat: (
      state,
      action: PayloadAction<{
        itemId: number
        saveFormat: string
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      targetItem.saveFormat = action.payload.saveFormat
    },
    setSaveFileName: (
      state,
      action: PayloadAction<{
        itemId: number
        saveFileName: string
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      targetItem.saveFileName = action.payload.saveFileName
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
      if (isImageItem(targetItem)) {
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
      if (isTimeSeriesItem(targetItem)) {
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
      if (isHeatMapItem(targetItem)) {
        targetItem.filePath = filePath
        targetItem.nodeId = nodeId
      }
    },
    resetImageActiveIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        startIndex?: number
        endIndex?: number
      }>,
    ) => {
      resetImageActiveIndexFn(state, action.payload)
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
    setImageActiveIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        activeIndex: number
      }>,
    ) => {
      const { itemId, activeIndex } = action.payload
      const targetItem = state.items[itemId]
      if (isImageItem(targetItem)) {
        targetItem.activeIndex = activeIndex
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
      }
    },
    setImageItemAlpha: (
      state,
      action: PayloadAction<{
        itemId: number
        alpha: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.alpha = action.payload.alpha
      }
    },
    setImageItemRoiAlpha: (
      state,
      action: PayloadAction<{
        itemId: number
        roiAlpha: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.roiAlpha = action.payload.roiAlpha
      }
    },
    setImageItemDuration: (
      state,
      action: PayloadAction<{
        itemId: number
        duration: number
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isImageItem(targetItem)) {
        targetItem.duration = action.payload.duration
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
    setTimeSeriesItemDrawOrderList: (
      state,
      action: PayloadAction<{
        itemId: number
        drawOrderList: string[]
      }>,
    ) => {
      const { itemId, drawOrderList } = action.payload
      const targetItem = state.items[itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.drawOrderList = drawOrderList
      }
    },
    setTimeSeriesItemMaxIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        maxIndex: number
      }>,
    ) => {
      const { itemId, maxIndex } = action.payload
      const targetItem = state.items[itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.maxIndex = maxIndex
      }
    },
    setTimeSeriesRefImageItemId: (
      state,
      action: PayloadAction<{
        itemId: number
        refImageItemId: number | null
      }>,
    ) => {
      const { itemId, refImageItemId } = action.payload
      const targetItem = state.items[itemId]
      if (isTimeSeriesItem(targetItem)) {
        targetItem.refImageItemId = refImageItemId ?? null
        targetItem.drawOrderList = []
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
    setCsvItemSetHeader: (
      state,
      action: PayloadAction<{
        itemId: number
        setHeader: number | null
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isCsvItem(targetItem)) {
        targetItem.setHeader = action.payload.setHeader
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
        xIndex: string
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
        yIndex: string
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isScatterItem(targetItem)) {
        targetItem.yIndex = action.payload.yIndex
      }
    },
    setBarItemIndex: (
      state,
      action: PayloadAction<{
        itemId: number
        index: string
      }>,
    ) => {
      const targetItem = state.items[action.payload.itemId]
      if (isBarItem(targetItem)) {
        targetItem.index = action.payload.index
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(deleteDisplayItem, (state, action) => {
        const itemId = action.payload.itemId

        if (isImageItem(state.items[itemId])) {
          Object.values(state.items).forEach((item) => {
            if (isTimeSeriesItem(item) && item.refImageItemId === itemId) {
              item.refImageItemId = null
            }
          })
        }

        delete state.items[itemId]
        if (itemId === state.selectedItemId) {
          state.selectedItemId = null
        }
        state.layout.forEach((row, i) => {
          const index = row.indexOf(itemId)
          if (index !== -1) {
            row.splice(index, 1)
          }
          if (row.length === 0) {
            state.layout.splice(i, 1)
          }
        })
      })
      .addCase(setNewDisplayDataPath, (state, action) => {
        const { itemId, filePath, nodeId, dataType } = action.payload
        const targetItem = state.items[itemId]
        if (isDisplayDataItem(targetItem)) {
          if (dataType != null) {
            state.items[itemId] = {
              ...getDisplayDataItemInitialValue(dataType),
              width: targetItem.width,
              height: targetItem.height,
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
        resetImageActiveIndexFn(state, { itemId })
      })
      .addCase(setImageItemClikedDataId.fulfilled, (state, action) => {
        const { itemId: imageItemId, clickedDataId } = action.meta.arg
        const targetItem = state.items[imageItemId]
        if (isImageItem(targetItem)) {
          targetItem.clickedDataId = clickedDataId
        }
        Object.values(state.items).forEach((item) => {
          if (isTimeSeriesItem(item)) {
            if (
              item.refImageItemId != null &&
              imageItemId === item.refImageItemId &&
              !item.drawOrderList.includes(clickedDataId)
            ) {
              item.drawOrderList.push(clickedDataId)
            }
          }
        })
      })
      .addCase(selectingImageArea.fulfilled, (state, action) => {
        const { itemId: imageItemId } = action.meta.arg
        const selectedZList = action.payload
        Object.values(state.items).forEach((item) => {
          if (isTimeSeriesItem(item)) {
            if (
              item.refImageItemId != null &&
              imageItemId === item.refImageItemId
            ) {
              item.drawOrderList = selectedZList
            }
          }
        })
      })
  },
})

function getMaxItemId(state: VisualaizeItem) {
  const idList = Object.keys(state.items).map((key) => Number(key))
  const maxId =
    idList.length > 0 ? idList.reduce((a, b) => Math.max(a, b)) : null
  return maxId
}

function addInitialItemFn(state: VisualaizeItem) {
  const maxId = getMaxItemId(state)
  const nextId = maxId != null ? maxId + 1 : 0
  state.items[nextId] = getDisplayDataItemInitialValue(DATA_TYPE_SET.IMAGE)
  state.selectedItemId = nextId
  return nextId
}

function resetImageActiveIndexFn(
  state: VisualaizeItem,
  args: {
    itemId: number
  },
) {
  const { itemId } = args
  const targetItem = state.items[itemId]
  if (isImageItem(targetItem)) {
    targetItem.activeIndex = 0
  }
}

export const {
  pushInitialItemToNewRow,
  insertInitialItemToNextColumn,
  addItemForWorkflowDialog,
  deleteAllItemForWorkflowDialog,
  setItemSize,
  selectItem,
  setFilePath,
  setSaveFormat,
  setSaveFileName,
  setHeatMapItemFilePath,
  setImageItemFilePath,
  setTimeSeriesItemFilePath,
  setRoiItemFilePath,
  resetImageActiveIndex,
  incrementImageActiveIndex,
  decrementImageActiveIndex,
  setImageActiveIndex,
  setImageItemShowticklabels,
  setImageItemZsmooth,
  setImageItemShowLine,
  setImageItemShowGrid,
  setImageItemShowScale,
  setImageItemColors,
  setImageItemStartIndex,
  setImageItemEndIndex,
  setImageItemAlpha,
  setImageItemRoiAlpha,
  setImageItemDuration,
  setTimeSeriesItemOffset,
  setTimeSeriesItemSpan,
  setTimeSeriesItemShowGrid,
  setTimeSeriesItemShowLine,
  setTimeSeriesItemShowTickLabels,
  setTimeSeriesItemZeroLine,
  setTimeSeriesItemXrangeLeft,
  setTimeSeriesItemXrangeRight,
  setTimeSeriesItemDrawOrderList,
  setTimeSeriesItemMaxIndex,
  setTimeSeriesRefImageItemId,
  setHeatMapItemShowScale,
  setHeatMapItemColors,
  setCsvItemTranspose,
  setCsvItemSetHeader,
  setCsvItemSetIndex,
  setScatterItemXIndex,
  setScatterItemYIndex,
  setBarItemIndex,
} = visualaizeItemSlice.actions

export default visualaizeItemSlice.reducer
