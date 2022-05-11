import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import { Node } from 'react-flow-renderer'
import { ParamDTO } from 'utils/param/ParamType'
import { FILE_TYPE } from '../InputNode/InputNodeType'
import { FLOW_ELEMENT_SLICE_NAME, NodeData } from './FlowElementType'
import { getAlgoParamsApi } from 'api/params/Params'

export const addAlgorithmNode = createAsyncThunk<
  ParamDTO,
  {
    node: Omit<Node<NodeData>, 'position'>
    functionPath: string
    name: string
  }
>(`${FLOW_ELEMENT_SLICE_NAME}/addAlgorithmNode`, async ({ name }, thunkAPI) => {
  try {
    const response = await getAlgoParamsApi(name)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const addInputNode = createAction<{
  node: Omit<Node<NodeData>, 'position'>
  fileType: FILE_TYPE
}>(`${FLOW_ELEMENT_SLICE_NAME}/addInputNode`)
