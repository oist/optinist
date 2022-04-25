import { RootState } from 'store/store'
import { selectRoiData } from '../DisplayData/DisplayDataSelectors'
import { DATA_TYPE } from '../DisplayData/DisplayDataType'
import {
  isDisplayDataItem,
  isImageItem,
  isTimeSeriesItem,
  isHeatMapItem,
  isCsvItem,
  isScatterItem,
} from './VisualizeItemUtils'

export const selectSelectedVisualizeItemId = (state: RootState) =>
  state.visualaizeItem.selectedItemId

export const selectVisualizeImageItemIdList = (state: RootState) =>
  Object.keys(state.visualaizeItem.items)
    .map(Number)
    .filter((itemId) => {
      const item = selectVisualizeItems(state)[itemId]
      return isImageItem(item) && !item.isWorkflowDialog
    })

export const selectVisualizeItems = (state: RootState) =>
  state.visualaizeItem.items

export const selectVisualizeItemLayout = (state: RootState) =>
  state.visualaizeItem.layout

export const selectVisualizeItemIsWorkflowDialog =
  (itemId: number) => (state: RootState) => {
    return selectVisualizeItems(state)[itemId].isWorkflowDialog
  }

export const selectVisualizeItemIdForWorkflowDialog =
  (nodeId: string, filePath: string, dataType: DATA_TYPE) =>
  (state: RootState) => {
    const items = selectVisualizeItems(state)
    let targetItemId: null | number = null
    for (const [itemId, value] of Object.entries(items)) {
      if (
        value.nodeId === nodeId &&
        value.filePath === filePath &&
        value.dataType === dataType &&
        value.isWorkflowDialog
      ) {
        targetItemId = Number(itemId)
      }
    }
    return targetItemId
  }

export const selectVisualizeItemWidth =
  (itemId: number) => (state: RootState) => {
    return selectVisualizeItems(state)[itemId].width
  }

export const selectVisualizeItemHeight =
  (itemId: number) => (state: RootState) => {
    return selectVisualizeItems(state)[itemId].height
  }

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

export const selectImageItemFilePath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDisplayDataItem(item)) {
      return item.filePath
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemFilePath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDisplayDataItem(item)) {
      return item.filePath
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectRoiItemNodeId = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isImageItem(item)) {
    return item.roiItem?.nodeId ?? null
  } else {
    throw new Error('invalid VisualaizeItemType')
  }
}

export const selectRoiItemFilePath = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isImageItem(item)) {
    return item.roiItem?.filePath ?? null
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

export const selectImageItemStartIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.startIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemEndIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.endIndex
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

export const selectImageItemRoiAlpha =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.roiAlpha
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemDuration =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.duration
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemSaveFilename =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.saveFileName
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemSaveFormat =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.saveFormat
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

export const selectTimeSeriesItemDisplayNumbers =
  (itemId: number, refImageItemId?: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.displayNumbers
    }
    throw new Error('invalid VisualaizeItemType')
  }

export const selectTimeSeriesItemRefImageItemId =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.refImageItemId
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemCheckedList =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.checkedList
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectRoiItemIndex =
  (itemId: number, roiFilePath: string | null) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      const maxIdx = item.maxIndex
      if (maxIdx !== 0) {
        return maxIdx
      }
    }

    if (roiFilePath !== null) {
      return selectRoiItemMaxNumber(roiFilePath)(state)
    } else {
      return 0
    }
  }

export const selectRoiItemMaxNumber =
  (roiFilePath: string) => (state: RootState) => {
    const roiData = selectRoiData(roiFilePath)(state)
    if (roiData.length !== 0) {
      return Math.max(...roiData.map((arr) => Math.max(...arr)))
    } else {
      return 0
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

export const selectCsvItemTranspose =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isCsvItem(item)) {
      return item.transpose
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectCsvItemSetHeader =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isCsvItem(item)) {
      return item.setHeader
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectCsvItemSetIndex = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isCsvItem(item)) {
    return item.setIndex
  } else {
    throw new Error('invalid VisualaizeItemType')
  }
}

export const selectScatterItemXIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isScatterItem(item)) {
      return item.xIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectScatterItemYIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isScatterItem(item)) {
      return item.yIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectDisplayDataIsSingle =
  (itemId: number) => (state: RootState) => {
    const items = selectVisualizeItems(state)
    const item = items[itemId]
    const targetFilePath = item.filePath
    return (
      Object.values(items).filter((value) => value.filePath === targetFilePath)
        .length === 1
    )
  }
