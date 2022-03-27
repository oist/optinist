import { createAction, createAsyncThunk } from '@reduxjs/toolkit'
import { Node } from 'react-flow-renderer'
import axios from 'axios'
import { BASE_URL } from 'const/API'
import { ParamDTO } from 'store/utils/param/ParamType'
import { FILE_TYPE } from '../InputNode/InputNodeType'
import { FLOW_ELEMENT_SLICE_NAME, NodeData } from './FlowElementType'

export const addAlgorithmNode = createAsyncThunk<
  ParamDTO,
  {
    node: Omit<Node<NodeData>, 'position'>
    functionPath: string
    name: string
  }
>(`${FLOW_ELEMENT_SLICE_NAME}/addAlgorithmNode`, async ({ name }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/params/${name}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const addInputNode = createAction<{
  node: Omit<Node<NodeData>, 'position'>
  fileType: FILE_TYPE
}>(`${FLOW_ELEMENT_SLICE_NAME}/addInputNode`)
