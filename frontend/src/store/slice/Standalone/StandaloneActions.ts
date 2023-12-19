import { createAsyncThunk } from "@reduxjs/toolkit"

import { getModeStandaloneApi } from "api/modeStandalone"

export const getModeStandalone = createAsyncThunk<boolean>(
  "/getModeStandalone",
  async () => {
    const response = await getModeStandaloneApi()
    return response
  },
)
