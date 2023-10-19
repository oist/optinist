import { createAsyncThunk, createAction } from "@reduxjs/toolkit"


import { getTimeSeriesDataById } from "store/slice/DisplayData/DisplayDataActions"
import { selectRoiData } from "store/slice/DisplayData/DisplayDataSelectors"
import { DATA_TYPE } from "store/slice/DisplayData/DisplayDataType"
import { selectVisualizeItems } from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { VISUALIZE_ITEM_SLICE_NAME } from "store/slice/VisualizeItem/VisualizeItemType"
import { isImageItem, isTimeSeriesItem } from "store/slice/VisualizeItem/VisualizeItemUtils"
import { ThunkApiConfig } from "store/store"

export const setImageItemClikedDataId = createAsyncThunk<
  void,
  { itemId: number; clickedDataId: string },
  ThunkApiConfig
>(
  `${VISUALIZE_ITEM_SLICE_NAME}/setImageItemClikedDataId`,
  ({ itemId, clickedDataId }, thunkAPI) => {
    const items = selectVisualizeItems(thunkAPI.getState())
    Object.values(items).forEach((item) => {
      if (
        isTimeSeriesItem(item) &&
        item.filePath != null &&
        item.refImageItemId === itemId &&
        !item.drawOrderList.includes(clickedDataId)
      ) {
        thunkAPI.dispatch(
          getTimeSeriesDataById({ path: item.filePath, index: clickedDataId }),
        )
      }
    })
  },
)

export const selectingImageArea = createAsyncThunk<
  string[],
  {
    itemId: number
    range: {
      x: number[]
      y: number[]
    }
  },
  ThunkApiConfig
>(
  `${VISUALIZE_ITEM_SLICE_NAME}/selectingImageArea`,
  ({ itemId, range }, thunkAPI) => {
    const { x, y } = range
    const [x1, x2] = x.map(Math.round)
    const [y1, y2] = y.map(Math.round)
    const selectedZList: string[] = []
    const items = selectVisualizeItems(thunkAPI.getState())
    const imageItem = items[itemId]
    if (isImageItem(imageItem) && imageItem.roiItem != null) {
      const roiFilePath = imageItem.roiItem.filePath
      if (roiFilePath != null) {
        const roiData = selectRoiData(roiFilePath)(thunkAPI.getState())
        for (let x = x1; x <= x2; x++) {
          for (let y = y1; y <= y2; y++) {
            const z = roiData[y][x]
            if (z != null) {
              const zNum = String(z) // indexとidのずれを回避
              if (!selectedZList.includes(zNum)) {
                selectedZList.push(zNum)
              }
            }
          }
        }
        Object.values(items).forEach((item) => {
          if (
            isTimeSeriesItem(item) &&
            item.filePath != null &&
            item.refImageItemId === itemId
          ) {
            const path = item.filePath
            selectedZList.forEach((selectedZ) => {
              thunkAPI.dispatch(
                getTimeSeriesDataById({
                  path,
                  index: String(selectedZ),
                }),
              )
            })
          }
        })
      }
    }
    return selectedZList
  },
)

type DeleteOption =
  | {
      deleteData: true
      dataType: DATA_TYPE
      filePath: string
    }
  | { deleteData: false }

type DeleteDisplayItemArgs = {
  itemId: number
} & DeleteOption

export const deleteDisplayItem = createAction<DeleteDisplayItemArgs>(
  `${VISUALIZE_ITEM_SLICE_NAME}/deleteDisplayItem`,
)

type SetNewDisplayDataPathArgs = {
  itemId: number
  filePath: string
  nodeId: string | null
  dataType?: DATA_TYPE
} & (
  | {
      deleteData: true
      prevFilePath: string
      prevDataType: DATA_TYPE
    }
  | { deleteData: false }
)

export const setNewDisplayDataPath = createAction<SetNewDisplayDataPathArgs>(
  `${VISUALIZE_ITEM_SLICE_NAME}/setNewDisplayDataPath`,
)
