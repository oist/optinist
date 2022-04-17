import { createAsyncThunk, createAction } from '@reduxjs/toolkit'
import { ThunkApiConfig } from 'store/store'
import { getTimeSeriesDataById } from '../DisplayData/DisplayDataActions'
import { DATA_TYPE } from '../DisplayData/DisplayDataType'
import { selectVisualizeItems } from './VisualizeItemSelectors'
import { VISUALIZE_ITEM_SLICE_NAME } from './VisualizeItemType'
import { isTimeSeriesItem } from './VisualizeItemUtils'

export const setImageItemClikedDataId = createAsyncThunk<
  void,
  { itemId: number; clickedDataId: number },
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
        !item.displayNumbers.includes(clickedDataId)
      ) {
        thunkAPI.dispatch(
          getTimeSeriesDataById({ path: item.filePath, index: clickedDataId }),
        )
      }
    })
    return
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
