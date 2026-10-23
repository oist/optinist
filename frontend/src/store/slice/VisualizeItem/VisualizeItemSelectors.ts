import { createSelector } from "@reduxjs/toolkit"

import {
  selectRoiData,
  selectRoiUniqueList,
  selectTimeSeriesData,
  selectTimeSeriesDataIsInitialized,
} from "store/slice/DisplayData/DisplayDataSelectors"
import { DATA_TYPE } from "store/slice/DisplayData/DisplayDataType"
import {
  isDisplayDataItem,
  isImageItem,
  isTimeSeriesItem,
  isHeatMapItem,
  isCsvItem,
  isScatterItem,
  isBarItem,
  isHistogramItem,
  isLineItem,
  isPolarItem,
  isRoiItem,
} from "store/slice/VisualizeItem/VisualizeItemUtils"
import { RootState } from "store/store"

export const selectSelectedVisualizeItemId = (state: RootState) =>
  state.visualaizeItem.selectedItemId

export const selectVisualizeImageItemIdList = (state: RootState) =>
  Object.keys(state.visualaizeItem.items)
    .map(Number)
    .filter((itemId) => {
      const item = selectVisualizeItemById(itemId)(state)
      return isImageItem(item) && !item.isWorkflowDialog
    })

export const selectVisualizeImageAndRoiItemIdList = (state: RootState) =>
  Object.keys(state.visualaizeItem.items)
    .map(Number)
    .filter((itemId) => {
      const item = selectVisualizeItemById(itemId)(state)
      return !item.isWorkflowDialog && (isImageItem(item) || isRoiItem(item))
    })

export const selectVisualizeItems = (state: RootState) =>
  state.visualaizeItem.items

export const selectVisualizeItemById = (itemId: number) => (state: RootState) =>
  state.visualaizeItem.items[itemId]

export const selectVisualizeItemLayout = (state: RootState) =>
  state.visualaizeItem.layout

export const selectVisualizeItemIsWorkflowDialog =
  (itemId: number) => (state: RootState) => {
    return selectVisualizeItemById(itemId)(state).isWorkflowDialog
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
    return selectVisualizeItemById(itemId)(state).width
  }

export const selectVisualizeItemHeight =
  (itemId: number) => (state: RootState) => {
    return selectVisualizeItemById(itemId)(state).height
  }

export const selectVisualizeItemType = (itemId: number) => (state: RootState) =>
  selectVisualizeItemById(itemId)(state).itemType

export const selectVisualizeDataType =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isDisplayDataItem(item)) {
      return item.dataType
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectVisualizeDataNodeId =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isDisplayDataItem(item)) {
      return item.nodeId
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectVisualizeDataFilePath =
  (itemId: number | null) => (state: RootState) => {
    if (!itemId && itemId !== 0) return ""
    const item = selectVisualizeItemById(itemId)(state)
    if (isDisplayDataItem(item)) {
      return item.filePath
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectVisualizeSaveFilename =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    return item.saveFileName
  }

export const selectVisualizeSaveFormat =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    return item.saveFormat
  }

export const selectImageItemFilePath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isDisplayDataItem(item)) {
      return item.filePath
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemFilePath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isDisplayDataItem(item)) {
      return item.filePath
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectRoiItemNodeId = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItemById(itemId)(state)
  if (isImageItem(item)) {
    return item.roiItem?.nodeId ?? null
  } else {
    throw new Error("invalid VisualaizeItemType")
  }
}

export const selectRoiItemFilePath = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItemById(itemId)(state)
  if (isImageItem(item)) {
    return item.roiItem?.filePath ?? null
  } else {
    throw new Error("invalid VisualaizeItemType")
  }
}

export const selectRoiItemOutputKeys =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.roiItem?.outputKey ?? null
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemShowticklabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.showticklabels
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemZsmooth =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.zsmooth
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemStartIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.startIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemEndIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.endIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemShowLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.showline
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemShowGrid =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.showgrid
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemShowScale =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.showscale
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemColors = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItemById(itemId)(state)
  if (isImageItem(item)) {
    return item.colors
  } else {
    throw new Error("invalid VisualaizeItemType")
  }
}

export const selectImageItemActiveIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.activeIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemAlpha = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItemById(itemId)(state)
  if (isImageItem(item)) {
    return item.alpha
  } else {
    throw new Error("invalid VisualaizeItemType")
  }
}

export const selectImageItemRoiAlpha =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.roiAlpha
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemRoiFilePath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      if (item.roiItem) {
        return item.roiItem.filePath
      } else {
        return null
      }
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemDuration =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.duration
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemOffset =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.stdBool
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemSpan =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.span
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemShowGrid =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.showgrid
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemShowLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.showline
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemShowTickLabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.showticklabels
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemZeroLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.zeroline
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemXrange =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.xrange
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemDrawOrderList =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.drawOrderList
    }
    throw new Error("invalid VisualaizeItemType")
  }

export const selectTimeSeriesItemRefImageItemId =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.refImageItemId
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemRefRoiUniqueList =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      if (item.refImageItemId != null) {
        const imageItem = selectVisualizeItems(state)[item.refImageItemId]
        if (isImageItem(imageItem) && imageItem.roiItem?.filePath) {
          return selectRoiUniqueList(imageItem.roiItem?.filePath)(state)
        }
        if (isRoiItem(imageItem) && imageItem.filePath) {
          return selectRoiUniqueList(imageItem.filePath)(state)
        }
      }
      return null
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectTimeSeriesItemKeys = (itemId: number) =>
  createSelector(
    [
      (state: RootState) => selectVisualizeItemById(itemId)(state),
      (state: RootState) => selectTimeSeriesItemFilePath(itemId)(state),
      (state: RootState) => {
        const path = selectTimeSeriesItemFilePath(itemId)(state)
        return path != null && selectTimeSeriesDataIsInitialized(path)(state)
      },
      (state: RootState) => {
        const path = selectTimeSeriesItemFilePath(itemId)(state)
        return path != null ? selectTimeSeriesData(path)(state) : {}
      },
      (state: RootState) => selectTimeSeriesItemRefRoiUniqueList(itemId)(state),
    ],
    (item, path, isInitialized, timeSeriesData, roiUniqueList) => {
      if (!isTimeSeriesItem(item)) {
        throw new Error("invalid VisualaizeItemType")
      }
      if (path != null && isInitialized) {
        const dataKeys = Object.keys(timeSeriesData)
        if (roiUniqueList != null) {
          return dataKeys.filter((key) => roiUniqueList.includes(key))
        } else {
          return dataKeys
        }
      } else {
        return []
      }
    },
  )

export const selectRoiItemIndex =
  (itemId: number, roiFilePath: string | null) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
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

export const selectRoiItemMaxNumber = (roiFilePath: string) =>
  createSelector(
    [(state: RootState) => selectRoiData(roiFilePath)(state)],
    (roiData) => {
      if (roiData.length !== 0) {
        return Math.max(...roiData.map((arr) => Math.max(...arr)))
      } else {
        return 0
      }
    },
  )

export const selectHeatMapItemShowScale =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isHeatMapItem(item)) {
      return item.showscale
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectHeatMapItemColors =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isHeatMapItem(item)) {
      return item.colors
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectCsvItemTranspose =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isCsvItem(item)) {
      return item.transpose
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectCsvItemSetHeader =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isCsvItem(item)) {
      return item.setHeader
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectCsvItemSetIndex = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItemById(itemId)(state)
  if (isCsvItem(item)) {
    return item.setIndex
  } else {
    throw new Error("invalid VisualaizeItemType")
  }
}

export const selectMatlabItemTranspose =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isCsvItem(item)) {
      return item.transpose
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectMatlabItemSetHeader =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isCsvItem(item)) {
      return item.setHeader
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectMatlabItemSetIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isCsvItem(item)) {
      return item.setIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectBarItemIndex = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItemById(itemId)(state)
  if (isBarItem(item)) {
    return item.index
  } else {
    throw new Error("invalid VisualaizeItemType")
  }
}

export const selectScatterItemXIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isScatterItem(item)) {
      return item.xIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectScatterItemYIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isScatterItem(item)) {
      return item.yIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
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

export const selectHistogramItemBins =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isHistogramItem(item)) {
      return item.bins
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectLineItemSelectedIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isLineItem(item)) {
      return item.selectedIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectPolarItemSelectedIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isPolarItem(item)) {
      return item.selectedIndex
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectImageItemRangeUnit =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isTimeSeriesItem(item)) {
      return item.rangeUnit
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }

export const selectClickedRoi = (itemId: number) => (state: RootState) =>
  state.visualaizeItem.clickedRois[itemId] || null

export const selectImageItemShowRoiLabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItemById(itemId)(state)
    if (isImageItem(item)) {
      return item.showRoiLabels
    } else {
      throw new Error("invalid VisualaizeItemType")
    }
  }
