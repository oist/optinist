import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { ELEMENT_SLICE_NAME } from './ElementType'
import { NODE_DATA_TYPE } from 'const/NodeData'
import { ThunkApiConfig } from 'redux/store'
import { flowElementsSelector } from './ElementSelector'
import { isInputNodeData, isAlgoNodeData } from 'utils/ElementUtils'
import { OutputData } from '../Algorithm/AlgorithmType'

export const clickNode = createAction<{ id: string; type: NODE_DATA_TYPE }>(
  `${ELEMENT_SLICE_NAME}/clickNode`,
)

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
