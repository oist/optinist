import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { ELEMENT_SLICE_NAME, Param } from './ElementType'
import { NODE_DATA_TYPE } from 'const/NodeData'
import { ThunkApiConfig } from 'redux/store'
import { flowElementsSelector } from './ElementSelector'
import { isInputNodeData, isAlgoNodeData } from 'utils/ElementUtils'
import { OutputData } from '../Output/OutputType'

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

export const runPipeline = createAsyncThunk<
  { message: string; data: OutputData[] },
  void,
  ThunkApiConfig
>(`${ELEMENT_SLICE_NAME}/runPipeline`, async (_, thunkAPI) => {
  const nodeDataListForRun = flowElementsSelector(thunkAPI.getState())
    .filter((element) => isInputNodeData(element) || isAlgoNodeData(element))
    .map((element) => element.data)
  try {
    const response = await axios.post(
      'http://localhost:8000/run',
      nodeDataListForRun,
    )
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const stopPipeline = createAction<void>(
  `${ELEMENT_SLICE_NAME}/stopPipeline`,
)
