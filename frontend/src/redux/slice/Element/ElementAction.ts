import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { ELEMENT_SLICE_NAME } from './ElementType'
import { AlgoNodeData, NODE_DATA_TYPE } from 'const/NodeData'
import { ThunkApiConfig } from 'redux/store'
import { flowElementsSelector } from './ElementSelector'
import { isInputNodeData, isAlgoNodeData } from 'utils/ElementUtils'
import { algoParamByIdSelector } from '../Algorithm/AlgorithmSelector'
import { BASE_URL } from 'const/API'

export const clickNode = createAction<{ id: string; type: NODE_DATA_TYPE }>(
  `${ELEMENT_SLICE_NAME}/clickNode`,
)

type OutputPathsDTO = {
  [algoName: string]: {
    [key: string]: {
      path: string
      type: string
      max_index?: number
    }
  }
}

export const runPipeline = createAsyncThunk<
  { message: string; outputPaths: OutputPathsDTO },
  void,
  ThunkApiConfig
>(`${ELEMENT_SLICE_NAME}/runPipeline`, async (_, thunkAPI) => {
  const pathErrorNodeList = flowElementsSelector(thunkAPI.getState()).filter(
    (element) => {
      if (isInputNodeData(element)) {
        if (!element.data?.path) {
          return true
        }
      }
      return false
    },
  )
  if (pathErrorNodeList.length > 0) {
    const message = `failed to read file path.`
    console.log(message, JSON.stringify(pathErrorNodeList))
    throw new Error(`${message}`)
  }
  const nodeDataListForRun = flowElementsSelector(thunkAPI.getState())
    .filter((element) => isInputNodeData(element) || isAlgoNodeData(element))
    .map((element) => {
      if (element.data && element.data.type === 'algo') {
        const param = algoParamByIdSelector(element.id)(thunkAPI.getState())
        const data: AlgoNodeData = {
          ...element.data,
          param,
        }
        return data
      } else {
        return element.data
      }
    })
  try {
    const response = await axios.post(`${BASE_URL}/api/run`, nodeDataListForRun)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const stopPipeline = createAction<void>(
  `${ELEMENT_SLICE_NAME}/stopPipeline`,
)
