import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import Select, { SelectChangeEvent } from '@mui/material/Select'
import FormControl from '@mui/material/FormControl'

import { arrayEqualityFn } from 'utils/EqualityUtils'
import { RootState } from 'store/store'
import {
  selectDisplayDataIsSingle,
  selectImageItemFilePath,
  selectRoiItemFilePath,
  selectRoiItemNodeId,
  selectSelectedVisualizeItemId,
  selectTimeSeriesItemRefImageItemId,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeImageItemIdList,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  selectItem,
  setItemSize,
  setRoiItemFilePath,
  setTimeSeriesRefImageItemId,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import { setNewDisplayDataPath } from 'store/slice/VisualizeItem/VisualizeItemActions'
import { useMouseDragHandler } from 'components/utils/MouseDragUtil'
import { DisplayDataItemLayoutMenuIcon } from './VisualizeItemLayoutMenuIcon'
import { DisplayDataItem } from './DisplayDataItem'
import { FilePathSelect } from './FilePathSelect'

export const VisualizeItem = React.memo<{ itemId: number }>(({ itemId }) => {
  const dispatch = useDispatch()
  const onClick = () => {
    dispatch(selectItem(itemId))
  }
  const isSelected = useSelector(
    (state: RootState) => selectSelectedVisualizeItemId(state) === itemId,
  )
  const { size, onMouseDownX, onMouseDownY, onMouseDownXY } =
    useItemDragResize(itemId)
  return (
    <Box sx={{ m: 1, display: 'flex', flexDirection: 'row' }}>
      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
        }}
      >
        <Paper
          variant="outlined"
          key={itemId}
          onClick={onClick}
          sx={{
            width: `${size.width}px`,
            height: `${size.height}px`,
            p: 1,
            borderColor: (theme) =>
              isSelected ? theme.palette.primary.light : undefined,
          }}
        >
          <ItemHeader itemId={itemId} />
          <DisplayDataItem itemId={itemId} />
        </Paper>
        <Box
          sx={{
            display: 'flex',
          }}
        >
          <Box
            sx={{
              flexGrow: 1,
              position: 'relative',
              top: '-2px',
              height: '4px',
              cursor: 'row-resize',
            }}
            onMouseDown={onMouseDownY}
          />
          <Box
            sx={{
              position: 'relative',
              top: '-2px',
              height: '4px',
              width: '12px',
              cursor: 'nwse-resize',
            }}
            onMouseDown={onMouseDownXY}
          />
        </Box>
      </Box>
      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
        }}
      >
        <Box
          sx={{
            flexGrow: 1,
            position: 'relative',
            left: '-2px',
            width: '4px',
            cursor: 'col-resize',
          }}
          onMouseDown={onMouseDownX}
        />
        <Box
          sx={{
            position: 'relative',
            height: '12px',
            width: '4px',
            left: '-2px',
            cursor: 'nwse-resize',
          }}
          onMouseDown={onMouseDownXY}
        />
      </Box>
    </Box>
  )
})

const ItemHeader = React.memo<{ itemId: number }>(({ itemId }) => {
  const itemDataType = useSelector(selectVisualizeDataType(itemId))
  return (
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
      {itemDataType === DATA_TYPE_SET.IMAGE && (
        <Box flexGrow={1}>
          <RoiSelect itemId={itemId} />
        </Box>
      )}
      <Box>
        <DisplayDataItemLayoutMenuIcon itemId={itemId} />
      </Box>
    </Box>
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

const RefImageItemIdSelect = React.memo<{
  itemId: number
}>(({ itemId }) => {
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

const MIN_WIDTH = 200
const MIN_HEIFHT = 150

function useItemDragResize(itemId: number) {
  const dispatch = useDispatch()
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const [movingSize, setMovingSize] = React.useState({ width, height })
  const onCommitSize = React.useCallback(
    (size: { width: number; height: number }) =>
      dispatch(setItemSize({ itemId, ...size })),
    [dispatch, itemId],
  )
  const onMouseDownX = useMouseDragHandler(
    (downEvent) => {
      let movingX = downEvent.screenX
      let movingWidth = movingSize.width
      return {
        onMouseMove: (moveEvent) => {
          const newWidth = movingWidth + (moveEvent.screenX - movingX)
          movingWidth = newWidth >= MIN_WIDTH ? newWidth : MIN_WIDTH
          setMovingSize((size) => ({ ...size, width: movingWidth }))
          movingX = moveEvent.screenX
        },
        onMouseUp: () => {
          onCommitSize({ ...movingSize, width: movingWidth })
        },
      }
    },
    [movingSize, onCommitSize],
  )
  const onMouseDownY = useMouseDragHandler(
    (downEvent) => {
      let movingY = downEvent.screenY
      let movingHeight = movingSize.height
      return {
        onMouseMove: (moveEvent) => {
          const newHeight = movingHeight + (moveEvent.screenY - movingY)
          movingHeight = newHeight >= MIN_HEIFHT ? newHeight : MIN_HEIFHT
          setMovingSize((size) => ({ ...size, height: movingHeight }))
          movingY = moveEvent.screenY
        },
        onMouseUp: () => {
          onCommitSize({ ...movingSize, height: movingHeight })
        },
      }
    },
    [movingSize, onCommitSize],
  )
  const onMouseDownXY = useMouseDragHandler(
    (downEvent) => {
      let movingX = downEvent.screenX
      let movingWidth = movingSize.width
      let movingY = downEvent.screenY
      let movingHeight = movingSize.height
      return {
        onMouseMove: (moveEvent) => {
          const newWidth = movingWidth + (moveEvent.screenX - movingX)
          movingWidth = newWidth >= MIN_WIDTH ? newWidth : MIN_WIDTH
          const newHeight = movingHeight + (moveEvent.screenY - movingY)
          movingHeight = newHeight >= MIN_HEIFHT ? newHeight : MIN_HEIFHT
          setMovingSize({ width: movingWidth, height: movingHeight })
          movingX = moveEvent.screenX
          movingY = moveEvent.screenY
        },
        onMouseUp: () => {
          onCommitSize({ width: movingWidth, height: movingHeight })
        },
      }
    },
    [movingSize, onCommitSize],
  )
  return {
    size: movingSize,
    onMouseDownX,
    onMouseDownY,
    onMouseDownXY,
  }
}

const RoiSelect = React.memo<{
  itemId: number
}>(({ itemId }) => {
  const dispatch = useDispatch()
  const roiItemNodeId = useSelector(selectRoiItemNodeId(itemId))
  const roiItemFilePath = useSelector(selectRoiItemFilePath(itemId))
  const onSelectRoiFilePath = (nodeId: string, filePath: string) => {
    dispatch(setRoiItemFilePath({ itemId, nodeId, filePath }))
  }
  return (
    <FilePathSelect
      selectedFilePath={roiItemFilePath}
      selectedNodeId={roiItemNodeId}
      onSelect={onSelectRoiFilePath}
      dataType={DATA_TYPE_SET.ROI}
      label={'Select Roi'}
    />
  )
})
