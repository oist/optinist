import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import TextField from '@mui/material/TextField'
import InputAdornment from '@mui/material/InputAdornment'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectSelectedVisualizeItemId,
  selectVisualizeItemLayout,
  selectVisualizeItemType,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { VisualizeItemAddButton } from './VisualizeItemAddButton'
import { VisualizeItemLayoutMenuIcon } from './VisualizeItemLayoutMenuIcon'
import { MultiPlotItem } from './MultiPlotItem'
import { DisplayDataItem } from './DisplayDataItem'
import {
  selectItem,
  setItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'

export const VisualizeItems: React.FC = () => {
  return (
    <>
      <FlexItemList />
    </>
  )
}

const FlexItemList: React.FC = () => {
  const layout = useSelector(selectVisualizeItemLayout, twoDimarrayEqualityFn)
  return (
    <Box display="flex" flexWrap="wrap" flexDirection="column" p={1} m={1}>
      {layout.map((row) => (
        <Box display="flex" flexDirection="row">
          {row.map((itemId) => (
            <Item itemId={itemId} key={itemId} />
          ))}
        </Box>
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

  const width = useSelector(selectVisualizeItemWidth(itemId))
  const [inputWidth, setInputWidth] = React.useState(width)
  const onBlurWidth = () => {
    const value = inputWidth >= 150 ? inputWidth : 150
    dispatch(
      setItemWidth({
        itemId,
        width: value,
      }),
    )
    setInputWidth(value)
  }
  const onChangeWidth = (event: React.ChangeEvent<HTMLInputElement>) => {
    const value = Number(event.target.value)
    setInputWidth(value)
  }
  return (
    <Paper
      variant="outlined"
      key={itemId}
      style={{
        width: `${width}px`,
        margin: theme.spacing(1),
        padding: theme.spacing(1),
        cursor: 'pointer',
        borderColor: isSelected ? theme.palette.primary.light : undefined,
      }}
      onClick={onSelect}
    >
      <Box display="flex" justifyContent="flex-end">
        <Box flexGrow={1}>
          <TextField
            type="number"
            size="small"
            label="width"
            sx={{ marginRight: 1, marginBottom: 1, width: '100px' }}
            InputProps={{
              endAdornment: <InputAdornment position="end">px</InputAdornment>,
            }}
            inputProps={{
              min: 150,
            }}
            value={inputWidth}
            onBlur={onBlurWidth}
            onChange={onChangeWidth}
          />
        </Box>
        <Box>
          <VisualizeItemLayoutMenuIcon itemId={itemId} />
        </Box>
      </Box>
      <ItemByType itemType={itemType} itemId={itemId} />
    </Paper>
  )
})

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
