import { createAsyncThunk } from "@reduxjs/toolkit"

import { getModeStandaloneApi } from "api/modeStandalone"

export const getModeStandalone = createAsyncThunk<boolean>(
  "/getModeStandalone",
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await getModeStandaloneApi()
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
