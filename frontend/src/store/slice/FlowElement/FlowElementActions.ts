import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import { Node } from 'react-flow-renderer'
import { ParamDTO } from 'utils/param/ParamType'
import { FILE_TYPE } from '../InputNode/InputNodeType'
import { FLOW_ELEMENT_SLICE_NAME, NodeData } from './FlowElementType'
import { getAlgoParamsApi } from 'api/params/Params'

export type AddingNodeData = Omit<Node<NodeData>, 'position'> &
  Partial<Pick<Node<NodeData>, 'position'>>

export const addAlgorithmNode = createAsyncThunk<
  ParamDTO,
  {
    node: AddingNodeData
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
  node: AddingNodeData
  fileType: FILE_TYPE
}>(`${FLOW_ELEMENT_SLICE_NAME}/addInputNode`)
