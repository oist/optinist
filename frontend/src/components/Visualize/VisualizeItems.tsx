import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@material-ui/core/styles'
import Box from '@material-ui/core/Box'
import Paper from '@material-ui/core/Paper'
import IconButton from '@material-ui/core/IconButton'
import MoreHorizIcon from '@material-ui/icons/MoreHoriz'

import { arrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectSelectedVisualizeItemId,
  selectVisualizeItemIdList,
  selectVisualizeItemType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { VisualizeItemAddButton } from './VisualizeItemAddButton'
import { VisualizeItemDeleteButton } from './VisualizeItemDeleteButton'
import { DefaultSetItem } from './DefaultSetItem'
import { DisplayDataItem } from './DisplayDataItem'
import { selectItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'

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
        // width: 300,
        // height: 300,
        margin: theme.spacing(1),
        padding: theme.spacing(1),
        borderColor: isSelected ? theme.palette.primary.light : undefined,
      }}
    >
      <Box display="flex">
        <Box flexGrow={1}>
          <IconButton onClick={onSelect}>
            <MoreHorizIcon />
          </IconButton>
        </Box>
        <Box>
          <VisualizeItemDeleteButton itemId={itemId} />
        </Box>
      </Box>
      <ItemByType itemId={itemId} />
    </Paper>
  )
})

const ItemByType = React.memo<{ itemId: number }>(({ itemId }) => {
  const itemType = useSelector(selectVisualizeItemType(itemId))
  switch (itemType) {
    case 'defaultSet':
      return <DefaultSetItem itemId={itemId} />
    case 'displayData':
      return <DisplayDataItem itemId={itemId} />
  }
})
