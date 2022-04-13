import React, { useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import Select, { SelectChangeEvent } from '@mui/material/Select'
import FormControl from '@mui/material/FormControl'

import { arrayEqualityFn, twoDimarrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectImageItemFilePath,
  selectSelectedVisualizeItemId,
  selectTimeSeriesItemRefImageItemId,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeImageItemIdList,
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
  setItemSize,
  setTimeSeriesRefImageItemId,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { FilePathSelect } from './FilePathSelect'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
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

  const [prevItem, setPrevItem] = React.useState<{
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
