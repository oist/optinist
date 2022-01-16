import {
  VisualaizeItemType,
  DefaultSetItem,
  DisplayDataItem,
  VISUALIZE_ITEM_TYPE_SET,
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
