import { Node } from "reactflow"

import { createAction, createAsyncThunk } from "@reduxjs/toolkit"

import { getAlgoParamsApi } from "api/params/Params"
import { FLOW_ELEMENT_SLICE_NAME, NodeData } from "store/slice/FlowElement/FlowElementType"
import { FILE_TYPE } from "store/slice/InputNode/InputNodeType"
import { ParamDTO } from "utils/param/ParamType"



export type AddingNodeData = Omit<Node<NodeData>, "position"> &
  Partial<Pick<Node<NodeData>, "position">>

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
