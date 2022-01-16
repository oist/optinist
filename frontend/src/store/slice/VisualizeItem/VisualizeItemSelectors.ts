import { RootState } from 'store/store'
import { selectNodeLabelById } from '../FlowElement/FlowElementSelectors'
import { isDisplayDataItem } from './VisualizeItemUtils'

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

// export const selectVisualizeDataNodeLabel =
//   (itemId: number) => (state: RootState) => {
//     const item = selectVisualizeItems(state)[itemId]
//     if (isDisplayDataItem(item)) {
//       return item.nodeId != null
//         ? selectNodeLabelById(item.nodeId)(state)
//         : null
//     } else {
//       throw new Error('failed to read dataType')
//     }
//   }
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
