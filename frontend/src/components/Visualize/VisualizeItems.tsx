import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'

import { arrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectSelectedVisualizeItemId,
  selectVisualizeItemIdList,
  selectVisualizeItemType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { VisualizeItemAddButton } from './VisualizeItemAddButton'
import { DisplayItemDeleteButton } from './DisplayItemDeleteButton'
import { MultiPlotItem } from './MultiPlotItem'
import { DisplayDataItem } from './DisplayDataItem'
import { selectItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { MultiPlotDeleteButton } from './MultiPlotDeleteButton'

export const VisualizeItems: React.FC = () => {
  return (
    <>
      <FlexItemList />
    </>
  )
}

const FlexItemList: React.FC = () => {
  const itemIdList = useSelector(selectVisualizeItemIdList, arrayEqualityFn)
  return (
    <Box display="flex" flexWrap="wrap" p={1} m={1}>
      {itemIdList.map((itemId) => (
        <Item itemId={itemId} key={itemId} />
      ))}
      <VisualizeItemAddButton />
    </Box>
  )
}

const Item = React.memo<{ itemId: number }>(({ itemId }) => {
  const itemType = useSelector(selectVisualizeItemType(itemId))

  const dispatch = useDispatch()
  const onSelect = () => {
    dispatch(selectItem(itemId))
  }
  const isSelected = useSelector(
    (state: RootState) => selectSelectedVisualizeItemId(state) === itemId,
  )
  const theme = useTheme()

  return (
    <Paper
      variant="outlined"
      key={itemId}
      style={{
        width: '100%',
        margin: theme.spacing(1),
        padding: theme.spacing(1),
        cursor: 'pointer',
        borderColor: isSelected ? theme.palette.primary.light : undefined,
      }}
      onClick={onSelect}
    >
      <Box display="flex" justifyContent="flex-end">
        <Box>
          <DeleteButton itemType={itemType} itemId={itemId} />
        </Box>
      </Box>
      <ItemByType itemType={itemType} itemId={itemId} />
    </Paper>
  )
})

const DeleteButton: React.FC<{
  itemType: string
  itemId: number
}> = ({ itemType, itemId }) => {
  switch (itemType) {
    case 'displayData':
      return <DisplayItemDeleteButton itemId={itemId} />
    case 'MultiPlot':
      return <MultiPlotDeleteButton itemId={itemId} />
    default:
      throw new Error('itemType Error')
  }
}

const ItemByType = React.memo<{
  itemType: string
  itemId: number
}>(({ itemType, itemId }) => {
  switch (itemType) {
    case 'MultiPlot':
      return <MultiPlotItem itemId={itemId} />
    case 'displayData':
      return <DisplayDataItem itemId={itemId} />
    default:
      throw new Error('itemType Error')
  }
})
