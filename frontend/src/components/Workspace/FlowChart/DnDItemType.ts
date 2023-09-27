export const DND_ITEM_TYPE_SET = {
  TREE_ITEM: 'TREE_ITEM',
} as const

export type DND_ITEM_TYPE =
  typeof DND_ITEM_TYPE_SET[keyof typeof DND_ITEM_TYPE_SET]

export type TreeItemDragObject = {}

export type TreeItemDropResult = {
  position?: { x: number; y: number }
}

export type TreeItemCollectedProps = {
  isDragging: boolean
}
