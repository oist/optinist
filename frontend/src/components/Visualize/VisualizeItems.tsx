import React, { useState, useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import TextField from '@mui/material/TextField'
import InputAdornment from '@mui/material/InputAdornment'
import { MenuItem, Select } from '@mui/material'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectImageItemFilePath,
  selectSelectedVisualizeItemId,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeIdList,
  selectVisualizeItemHeight,
  selectVisualizeItemLayout,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { VisualizeItemAddButton } from './VisualizeItemAddButton'
import { DisplayDataItemLayoutMenuIcon } from './VisualizeItemLayoutMenuIcon'
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

export const FlexItemList: React.FC = () => {
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

  const itemIdList = useSelector(selectVisualizeIdList)

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
          <>ID: {itemId}</>
          <FilePathSelectItem itemId={itemId} />
          <TextField
            type="number"
            size="small"
            label="width"
            sx={{ marginRight: 1, marginBottom: 1 }}
            InputProps={{
              endAdornment: <InputAdornment position="end">px</InputAdornment>,
            }}
            inputProps={{
              min: 150,
            }}
            style={{ width: 80, marginLeft: 10 }}
            value={inputWidth}
            onBlur={onBlurWidth}
            onChange={onChangeWidth}
          />
          <TextField
            type="number"
            label="height"
            sx={{ width: 80, marginRight: 1, marginBottom: 1 }}
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
        <Box flexGrow={1}>
          <Select label="smooth">
            {itemIdList.map((value) => (
              <MenuItem value={value}>{value}</MenuItem>
            ))}
          </Select>
        </Box>
        <Box>
          <DisplayDataItemLayoutMenuIcon itemId={itemId} />
        </Box>
      </Box>
      <DisplayDataItem itemId={itemId} />
    </Paper>
  )
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
