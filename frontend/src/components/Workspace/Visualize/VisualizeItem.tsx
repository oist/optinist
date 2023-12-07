import { memo, useCallback, useState, MouseEvent } from "react"
import { useDispatch, useSelector } from "react-redux"

import { Close, Numbers } from "@mui/icons-material"
import { Chip, IconButton } from "@mui/material"
import Box from "@mui/material/Box"
import FormControl from "@mui/material/FormControl"
import InputLabel from "@mui/material/InputLabel"
import MenuItem from "@mui/material/MenuItem"
import Paper from "@mui/material/Paper"
import Select, { SelectChangeEvent } from "@mui/material/Select"

import Loading from "components/common/Loading"
import { useMouseDragHandler } from "components/utils/MouseDragUtil"
import { DisplayDataItem } from "components/Workspace/Visualize/DisplayDataItem"
import { FilePathSelect } from "components/Workspace/Visualize/FilePathSelect"
import { selectLoadingVisualize } from "store/slice/DisplayData/DisplayDataSelectors"
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from "store/slice/DisplayData/DisplayDataType"
import {
  deleteDisplayItem,
  setNewDisplayDataPath,
} from "store/slice/VisualizeItem/VisualizeItemActions"
import {
  selectDisplayDataIsSingle,
  selectImageItemFilePath,
  selectRoiItemFilePath,
  selectRoiItemNodeId,
  selectSelectedVisualizeItemId,
  selectTimeSeriesItemRefImageItemId,
  selectVisualizeDataFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeImageItemIdList,
  selectVisualizeItemHeight,
  selectVisualizeItemWidth,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  selectItem,
  setItemSize,
  setRoiItemFilePath,
  setTimeSeriesRefImageItemId,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import { AppDispatch, RootState } from "store/store"
import { arrayEqualityFn } from "utils/EqualityUtils"

interface ItemIdProps {
  itemId: number
}

export const VisualizeItem = memo(function VisualizeItem({
  itemId,
}: ItemIdProps) {
  const dispatch = useDispatch()
  const onClick = () => {
    dispatch(selectItem(itemId))
  }
  const isSelected = useSelector(
    (state: RootState) => selectSelectedVisualizeItemId(state) === itemId,
  )
  const { size, onMouseDownX, onMouseDownY, onMouseDownXY } =
    useItemDragResize(itemId)
  const loading = useSelector(selectLoadingVisualize)
  return (
    <Box
      sx={{ m: 1, display: "flex", flexDirection: "row", position: "relative" }}
    >
      {loading ? <Loading position={"absolute"} /> : null}
      <Box
        sx={{
          display: "flex",
          flexDirection: "column",
        }}
      >
        <Paper
          variant="outlined"
          key={itemId}
          onClick={onClick}
          sx={{
            width: `${size.width}px`,
            minHeight: `${size.height}px`,
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
            display: "flex",
          }}
        >
          <Box
            sx={{
              flexGrow: 1,
              position: "relative",
              top: "-2px",
              height: "4px",
              cursor: "row-resize",
            }}
            onMouseDown={onMouseDownY}
          />
          <Box
            sx={{
              position: "relative",
              top: "-2px",
              height: "4px",
              width: "12px",
              cursor: "nwse-resize",
            }}
            onMouseDown={onMouseDownXY}
          />
        </Box>
      </Box>
      <Box
        sx={{
          display: "flex",
          flexDirection: "column",
        }}
      >
        <Box
          sx={{
            flexGrow: 1,
            position: "relative",
            left: "-2px",
            width: "4px",
            cursor: "col-resize",
          }}
          onMouseDown={onMouseDownX}
        />
        <Box
          sx={{
            position: "relative",
            height: "12px",
            width: "4px",
            left: "-2px",
            cursor: "nwse-resize",
          }}
          onMouseDown={onMouseDownXY}
        />
      </Box>
    </Box>
  )
})

const ItemHeader = memo(function ItemHeader({ itemId }: ItemIdProps) {
  const dataType = useSelector(selectVisualizeDataType(itemId))
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const isSingleData = useSelector(selectDisplayDataIsSingle(itemId))
  const dispatch = useDispatch<AppDispatch>()
  const handleClose = (e: MouseEvent) => {
    e.stopPropagation()
    dispatch(
      deleteDisplayItem(
        isSingleData && filePath != null && dataType != null
          ? { itemId, deleteData: true, filePath, dataType }
          : { itemId, deleteData: false },
      ),
    )
  }

  return (
    <Box display="flex" justifyContent="flex-end">
      <Box flexGrow={1} display="flex">
        <Chip
          icon={<Numbers />}
          size="small"
          label={itemId}
          color="primary"
          variant="outlined"
          sx={{ marginRight: 2 }}
        />
        <FilePathSelectItem itemId={itemId} />
      </Box>
      {dataType === DATA_TYPE_SET.TIME_SERIES && (
        <Box flexGrow={1}>
          <RefImageItemIdSelect itemId={itemId} />
        </Box>
      )}
      {dataType === DATA_TYPE_SET.IMAGE && (
        <Box flexGrow={1}>
          <RoiSelect itemId={itemId} />
        </Box>
      )}
      <Box>
        <IconButton onClick={handleClose}>
          <Close />
        </IconButton>
      </Box>
    </Box>
  )
})

const FilePathSelectItem = memo(function FilePathSelectItem({
  itemId,
}: ItemIdProps) {
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

const RefImageItemIdSelect = memo(function RefImageItemIdSelect({
  itemId,
}: ItemIdProps) {
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
        <MenuItem value={undefined}>{"None"}</MenuItem>
        {itemIdList.map((value) => (
          <MenuItem key={value} value={value}>
            {value}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
})

const MIN_WIDTH = 200
const MIN_HEIGHT = 150

function useItemDragResize(itemId: number) {
  const dispatch = useDispatch()
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const [movingSize, setMovingSize] = useState({ width, height })
  const onCommitSize = useCallback(
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
          movingHeight = newHeight >= MIN_HEIGHT ? newHeight : MIN_HEIGHT
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
          movingHeight = newHeight >= MIN_HEIGHT ? newHeight : MIN_HEIGHT
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

const RoiSelect = memo(function RoiSelect({ itemId }: ItemIdProps) {
  const dispatch = useDispatch()
  const roiItemNodeId = useSelector(selectRoiItemNodeId(itemId))
  const roiItemFilePath = useSelector(selectRoiItemFilePath(itemId))
  const onSelectRoiFilePath = (
    nodeId: string,
    filePath: string,
    dataType: string,
    outputKey?: string,
  ) => {
    dispatch(setRoiItemFilePath({ itemId, nodeId, filePath, outputKey }))
  }
  return (
    <FilePathSelect
      selectedFilePath={roiItemFilePath}
      selectedNodeId={roiItemNodeId}
      onSelect={onSelectRoiFilePath}
      dataType={DATA_TYPE_SET.ROI}
      label={"Select Roi"}
    />
  )
})
