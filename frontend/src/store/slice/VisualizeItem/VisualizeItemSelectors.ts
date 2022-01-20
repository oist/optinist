import { RootState } from 'store/store'
import {
  isDisplayDataItem,
  isImageItem,
  isTimeSeriesItem,
  isHeatMapItem,
} from './VisualizeItemUtils'

export const selectSelectedVisualizeItemId = (state: RootState) =>
  state.visualaizeItem.selectedItemId

const selectVisualizeItems = (state: RootState) => state.visualaizeItem.items

export const selectVisualizeItemIdList = (state: RootState) =>
  Object.keys(selectVisualizeItems(state)).map((key) => Number(key))

export const selectVisualizeItemType = (itemId: number) => (state: RootState) =>
  selectVisualizeItems(state)[itemId].itemType

export const selectVisualizeDataType =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDisplayDataItem(item)) {
      return item.dataType
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectVisualizeDataNodeId =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDisplayDataItem(item)) {
      return item.nodeId
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectVisualizeDataFilePath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDisplayDataItem(item)) {
      return item.filePath
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowticklabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showticklabels
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemZsmooth =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.zsmooth
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemMaxIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.maxIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showline
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowGrid =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showgrid
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowScale =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showscale
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemColors = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isImageItem(item)) {
    return item.colors
  } else {
    throw new Error('invalid VisualaizeItemType')
  }
}

export const selectImageItemActiveIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.activeIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemOffset =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.offset
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemSpan =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.span
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemShowGrid =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.showgrid
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemShowLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.showline
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemShowTickLabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.showticklabels
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemZeroLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.zeroline
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemXrange =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.xrange
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectHeatMapItemShowScale =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isHeatMapItem(item)) {
      return item.showscale
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectHeatMapItemColors =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isHeatMapItem(item)) {
      return item.colors
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }
