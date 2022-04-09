import { createAsyncThunk } from '@reduxjs/toolkit'
import { ThunkApiConfig } from 'store/store'
import { getTimeSeriesDataById } from '../DisplayData/DisplayDataActions'
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
