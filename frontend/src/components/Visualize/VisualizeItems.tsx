import React, { useState, useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import TextField from '@mui/material/TextField'
import InputAdornment from '@mui/material/InputAdornment'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectImageItemFilePath,
  selectSelectedVisualizeItemId,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeItemHeight,
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
  setDisplayDataPath,
  setItemHeight,
  setItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { FilePathSelect } from './FilePathSelect'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import { deleteDisplayItem } from 'store/slice/DisplayData/DisplayDataSlice'
import { VISUALIZE_ITEM_TYPE_SET } from 'store/slice/VisualizeItem/VisualizeItemType'

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
  const onClick = () => {
    dispatch(selectItem(itemId))
  }
  const isSelected = useSelector(
    (state: RootState) => selectSelectedVisualizeItemId(state) === itemId,
  )
  const theme = useTheme()

  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
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

  const [inputHeight, setInputHeight] = React.useState(height)
  const onBlurHeight = () => {
    const value = inputHeight >= 300 ? inputHeight : 300
    dispatch(
      setItemHeight({
        itemId,
        height: value,
      }),
    )
    setInputHeight(value)
  }
  const onChangeHeight = (event: React.ChangeEvent<HTMLInputElement>) => {
    const value = Number(event.target.value)
    setInputHeight(value)
  }

  return (
    <Paper
      variant="outlined"
      key={itemId}
      style={{
        width: `${width}px`,
        height: `${height}px`,
        margin: theme.spacing(1),
        padding: theme.spacing(1),
        cursor: 'pointer',
        borderColor: isSelected ? theme.palette.primary.light : undefined,
      }}
      onClick={onClick}
    >
      <Box display="flex" justifyContent="flex-end">
        <Box flexGrow={1}>
          {itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA && (
            <FilePathSelectItem itemId={itemId} />
          )}
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
            style={{ marginLeft: 10 }}
            value={inputWidth}
            onBlur={onBlurWidth}
            onChange={onChangeWidth}
          />
          <TextField
            type="number"
            size="small"
            label="height"
            sx={{ marginRight: 1, marginBottom: 1, width: '100px' }}
            InputProps={{
              endAdornment: <InputAdornment position="end">px</InputAdornment>,
            }}
            inputProps={{
              min: 150,
            }}
            value={inputHeight}
            onBlur={onBlurHeight}
            onChange={onChangeHeight}
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
    case VISUALIZE_ITEM_TYPE_SET.MULTI_PLOT:
      return <MultiPlotItem itemId={itemId} />
    case VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA:
      return <DisplayDataItem itemId={itemId} />
    default:
      throw new Error('itemType Error')
  }
})

const FilePathSelectItem = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dispatch = useDispatch()
  const dataType = useSelector(selectVisualizeDataType(itemId))
  const selectedNodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const selectedFilePath = useSelector(selectImageItemFilePath(itemId))

  const [prevItem, setPrevItem] = useState<{
    dataType: DATA_TYPE
    filePath: string | null
  }>({
    dataType: 'image',
    filePath: null,
  })

  useEffect(() => {
    setPrevItem({ dataType, filePath: selectedFilePath })
  }, [selectedFilePath, dataType])

  const onSelectFilePath = (
    nodeId: string,
    filePath: string,
    dataType: DATA_TYPE,
  ) => {
    dispatch(setDisplayDataPath({ itemId, nodeId, filePath, dataType }))
    dispatch(deleteDisplayItem(prevItem))
  }

  return (
    <FilePathSelect
      selectedNodeId={selectedNodeId}
      selectedFilePath={selectedFilePath}
      onSelect={onSelectFilePath}
    />
  )
})
