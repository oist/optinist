import { DATA_TYPE_SET } from '../DisplayData/DisplayDataType'
import {
  VisualaizeItemType,
  DefaultSetItem,
  DisplayDataItem,
  VISUALIZE_ITEM_TYPE_SET,
  ImageItem,
} from './VisualizeItemType'

export function isDefaultSetItem(
  item: VisualaizeItemType,
): item is DefaultSetItem {
  return item.itemType === VISUALIZE_ITEM_TYPE_SET.DEFAULT_SET
}

export function isDisplayDataItem(
  item: VisualaizeItemType,
): item is DisplayDataItem {
  return item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA
}

export function isImageItem(item: VisualaizeItemType): item is ImageItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.IMAGE
  )
}
