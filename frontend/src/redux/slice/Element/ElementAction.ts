import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { ELEMENT_SLICE_NAME, Param } from './ElementType'
import { NODE_DATA_TYPE } from 'const/NodeData'

export const clickNode = createAction<{ id: string; type: NODE_DATA_TYPE }>(
  `${ELEMENT_SLICE_NAME}/clickNode`,
)

export const getAlgoParams = createAsyncThunk<
  Param,
  { id: string; algoName: string }
>(`${ELEMENT_SLICE_NAME}/getAlgoParams`, async ({ algoName }, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await axios.get(`http://localhost:8000/params/${algoName}`)
    return response.data
  } catch (e) {
    return rejectWithValue(e)
  }
})
