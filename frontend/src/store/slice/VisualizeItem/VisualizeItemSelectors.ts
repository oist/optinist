import { RootState } from 'store/store'
import { VISUALIZE_ITEM_TYPE_SET } from './VisualizeItemType'
import {
  isDisplayDataItem,
  isImageItem,
  isTimeSeriesItem,
  isHeatMapItem,
  isRoiItem,
  isDefaultSetItem,
  isCsvItem,
  isScatterItem,
} from './VisualizeItemUtils'

export const selectSelectedVisualizeItemId = (state: RootState) =>
  state.visualaizeItem.selectedItemId

const selectVisualizeItems = (state: RootState) => state.visualaizeItem.items

export const selectVisualizeItemIdList = (state: RootState) =>
  Object.keys(selectVisualizeItems(state)).map((key) => Number(key))

export const selectVisualizeItemType = (itemId: number) => (state: RootState) =>
  selectVisualizeItems(state)[itemId].itemType

export const selectVisualizeItemTypeIsDefaultSet =
  (itemId: number) => (state: RootState) =>
    selectVisualizeItemType(itemId)(state) ===
    VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET

export const selectDefaultSetImageItem =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDefaultSetItem(item)) {
      return item.imageItem
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectDefaultSetImageItemNodeId =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetImageItem(itemId)(state).nodeId

export const selectDefaultSetImageItemFilePath =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetImageItem(itemId)(state).filePath

export const selectDefaultSetTimeSeriesItem =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDefaultSetItem(item)) {
      return item.timeSeriesItem
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectDefaultSetRoiItem = (itemId: number) => (state: RootState) =>
  selectDefaultSetImageItem(itemId)(state).roiItem

export const selectDefaultSetRoiItemNodeId =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetRoiItem(itemId)(state)?.nodeId ?? null

export const selectDefaultSetRoiItemFilePath =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetRoiItem(itemId)(state)?.filePath ?? null

export const selectDefaultSetTimeSeriesItemNodeId =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetTimeSeriesItem(itemId)(state).nodeId

export const selectDefaultSetTimeSeriesItemFilePath =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetTimeSeriesItem(itemId)(state).filePath

export const selectDefaultSetHeatMapItem =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDefaultSetItem(item)) {
      return item.heatMapItem
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectDefaultSetHeatMapItemNodeId =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetHeatMapItem(itemId)(state).nodeId

export const selectDefaultSetHeatMapItemFilePath =
  (itemId: number) => (state: RootState) =>
    selectDefaultSetHeatMapItem(itemId)(state).filePath

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
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.filePath
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectRoiItemNodeId = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isImageItem(item)) {
    return item.roiItem?.nodeId ?? null
  } else if (isDefaultSetItem(item)) {
    return item.imageItem.roiItem?.nodeId ?? null
  } else {
    throw new Error('invalid VisualaizeItemType')
  }
}

export const selectRoiItemFilePath = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isImageItem(item)) {
    return item.roiItem?.filePath ?? null
  } else if (isDefaultSetItem(item)) {
    return item.imageItem.roiItem?.filePath ?? null
  } else {
    throw new Error('invalid VisualaizeItemType')
  }
}

export const selectImageItemShowticklabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showticklabels
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.showticklabels
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemZsmooth =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.zsmooth
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.zsmooth
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemStartIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.startIndex
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.startIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemEndIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.endIndex
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.endIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showline
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.showline
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowGrid =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showgrid
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.showgrid
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemShowScale =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.showscale
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.showscale
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectImageItemColors = (itemId: number) => (state: RootState) => {
  const item = selectVisualizeItems(state)[itemId]
  if (isImageItem(item)) {
    return item.colors
  } else if (isDefaultSetItem(item)) {
    return item.imageItem.colors
  } else {
    throw new Error('invalid VisualaizeItemType')
  }
}

export const selectImageItemActiveIndex =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isImageItem(item)) {
      return item.activeIndex
    } else if (isDefaultSetItem(item)) {
      return item.imageItem.activeIndex
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemOffset =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.offset
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.offset
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemSpan =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.span
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.span
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemShowGrid =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.showgrid
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.showgrid
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemShowLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.showline
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.showline
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemShowTickLabels =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.showticklabels
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.showticklabels
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemZeroLine =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.zeroline
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.zeroline
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectTimeSeriesItemXrange =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.xrange
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.xrange
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

// export const selectTimeSeriesItemColorScale =
//   (itemId: number) => (state: RootState) => {
//     const item = selectVisualizeItems(state)[itemId]
//     if (isTimeSeriesItem(item)) {
//       return null
//     } else if (isDefaultSetItem(item)) {
//       const roiItem = item.imageItem.roiItem
//       if (roiItem !== null) {
//         return roiItem.colors
//       } else {
//         return null
//       }
//     } else {
//       throw new Error('invalid VisualaizeItemType')
//     }
//   }

export const selectTimeSeriesItemDisplayNumbers =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isTimeSeriesItem(item)) {
      return item.displayNumbers
    } else if (isDefaultSetItem(item)) {
      return item.timeSeriesItem.displayNumbers
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectHeatMapItemShowScale =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isHeatMapItem(item)) {
      return item.showscale
    } else if (isDefaultSetItem(item)) {
      return item.heatMapItem.showscale
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectHeatMapItemColors =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isHeatMapItem(item)) {
      return item.colors
    } else if (isDefaultSetItem(item)) {
      return item.heatMapItem.colors
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

// export const selectRoiItemColors = (itemId: number) => (state: RootState) => {
//   const item = selectVisualizeItems(state)[itemId]
//   if (isRoiItem(item)) {
//     return item.colors
//   } else if (isDefaultSetItem(item)) {
//     return item.imageItem.roiItem?.colors ?? []
//   } else if (isImageItem(item)) {
//     return item.roiItem?.colors ?? []
//   } else {
//     throw new Error('invalid VisualaizeItemType')
//   }
// }

export const selectCsvItemTranspose =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isCsvItem(item)) {
      return item.transpose
    } else {
      throw new Error('invalid VisualaizeItemType')
    }
  }

export const selectCsvItemSetColumn =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isCsvItem(item)) {
      return item.setColumn
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

export const selectDefaultSetTimeSeriesItemFilepath =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDefaultSetItem(item)) {
      const targetItem = selectDefaultSetTimeSeriesItem(itemId)(state)
      const timeSeriesFilePath = targetItem.filePath
      return timeSeriesFilePath
    }
    return null
  }

export const selectDefaultSetTimeSeriesItemDisplayNumbers =
  (itemId: number) => (state: RootState) => {
    const item = selectVisualizeItems(state)[itemId]
    if (isDefaultSetItem(item)) {
      const targetItem = selectDefaultSetTimeSeriesItem(itemId)(state)
      const displayNumbers = targetItem.displayNumbers
      return displayNumbers
    }
    return null
  }
