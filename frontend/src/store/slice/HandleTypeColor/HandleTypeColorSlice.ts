import { createSlice } from '@reduxjs/toolkit'
import * as MuiColors from '@mui/material/colors'

export type HandleTypeColor = {
  colorMap: { [type: string]: string }
  /**
   * HANDLE_COLOR_PRESET_MAPのkey
   */
  nextKey: number
}

const initialState: HandleTypeColor = {
  colorMap: {
    ImageData: MuiColors.red[500],
    IscellData: MuiColors.indigo[500],
    TimeSeriesData: MuiColors.yellow[500],
    Suite2pData: MuiColors.green[500],
    FluoData: MuiColors.orange[500],
    BehaviorData: MuiColors.yellow[500],
  },
  nextKey: 0,
}

const SLICE_NAME = 'handleTypeColor'

export const handleTypeColorSlice = createSlice({
  name: SLICE_NAME,
  initialState,
  reducers: {},
})

export default handleTypeColorSlice.reducer
