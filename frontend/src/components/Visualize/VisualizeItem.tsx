import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import Select, { SelectChangeEvent } from '@mui/material/Select'
import FormControl from '@mui/material/FormControl'

import { arrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectDisplayDataIsSingle,
  selectImageItemFilePath,
  selectSelectedVisualizeItemId,
  selectTimeSeriesItemRefImageItemId,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeImageItemIdList,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { DisplayDataItemLayoutMenuIcon } from './VisualizeItemLayoutMenuIcon'
import { DisplayDataItem } from './DisplayDataItem'
import {
  selectItem,
  setItemSize,
  setTimeSeriesRefImageItemId,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { FilePathSelect } from './FilePathSelect'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import { setNewDisplayDataPath } from 'store/slice/VisualizeItem/VisualizeItemActions'

export const VisualizeItem = React.memo<{ itemId: number }>(({ itemId }) => {
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

  const itemDataType = useSelector(selectVisualizeDataType(itemId))

  const [resizeTrigger, setResizeTrigger] = React.useState(false)
  const [resizeCoord, setResizeCoord] = React.useState<{
    x: number
    y: number
  }>({ x: 0, y: 0 })

  const onMouseDown = (event: React.MouseEvent<HTMLInputElement>) => {
    setResizeTrigger(true)
    setResizeCoord({ x: event.screenX, y: event.screenY })
  }

  const onMouseUp = () => {
    setResizeTrigger(false)
  }

  const onMouseLeave = () => {
    setResizeTrigger(false)
  }

  const onMouseMove = (event: React.MouseEvent<HTMLInputElement>) => {
    if (resizeTrigger) {
      const newWidth = width + (event.screenX - resizeCoord.x)
      const newHeight = height + (event.screenY - resizeCoord.y)
      dispatch(
        setItemSize({
          itemId,
          width: newWidth,
          height: newHeight,
        }),
      )
      setResizeCoord({ x: event.screenX, y: event.screenY })
    }
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
        cursor: resizeTrigger ? 'nwse-resize' : 'pointer',
        borderColor: isSelected ? theme.palette.primary.light : undefined,
      }}
      onClick={onClick}
      onMouseDown={onMouseDown}
      onMouseLeave={onMouseLeave}
      onMouseUp={onMouseUp}
      onMouseMove={onMouseMove}
    >
      <Box display="flex" justifyContent="flex-end">
        <Box flexGrow={1}>
          <>ID: {itemId}</>
          <FilePathSelectItem itemId={itemId} />
        </Box>
        {itemDataType === DATA_TYPE_SET.TIME_SERIES && (
          <Box flexGrow={1}>
            <RefImageItemIdSelect itemId={itemId} />
          </Box>
        )}
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

  const isSingleData = useSelector(selectDisplayDataIsSingle(itemId))
  const onSelectFilePath = (
    nodeId: string,
    newFilePath: string,
    newDataType: DATA_TYPE,
  ) => {
    const basePayload = {
      itemId,
      nodeId,
      filePath: newFilePath,
      dataType: newDataType,
    }
    dispatch(
      setNewDisplayDataPath(
        isSingleData && selectedFilePath != null
          ? {
              ...basePayload,
              deleteData: true,
              prevDataType: dataType,
              prevFilePath: selectedFilePath,
            }
          : {
              ...basePayload,
              deleteData: false,
            },
      ),
    )
  }

  return (
    <FilePathSelect
      selectedNodeId={selectedNodeId}
      selectedFilePath={selectedFilePath}
      onSelect={onSelectFilePath}
    />
  )
})

const RefImageItemIdSelect = React.memo<{ itemId: number }>(({ itemId }) => {
  const dispatch = useDispatch()
  const itemIdList = useSelector(
    selectVisualizeImageItemIdList,
    arrayEqualityFn,
  )
  const onChangeRefImageItemId = (event: SelectChangeEvent) => {
    const value = Number(event.target.value)
    dispatch(
      setTimeSeriesRefImageItemId({
        itemId,
        refImageItemId: isNaN(value) ? null : value,
      }),
    )
  }
  const selectedRefImageItemId = useSelector(
    selectTimeSeriesItemRefImageItemId(itemId),
  )
  return (
    <FormControl fullWidth variant="standard">
      <InputLabel>ref image</InputLabel>
      <Select
        value={String(selectedRefImageItemId)}
        onChange={onChangeRefImageItemId}
      >
        <MenuItem value={undefined}>{'None'}</MenuItem>
        {itemIdList.map((value) => (
          <MenuItem value={value}>{value}</MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})
