import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { ELEMENT_SLICE_NAME, NodeType, Param } from './ElementType'

export const clickNode = createAction<{ id: string; type: NodeType }>(
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
