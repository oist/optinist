import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { isInputNodePostData } from 'api/run/RunUtils'
import { INITIAL_IMAGE_ELEMENT_ID } from 'const/flowchart'
import {
  fetchExperiment,
  importExperimentByUid,
} from '../Experiments/ExperimentsActions'
import { uploadFile } from '../FileUploader/FileUploaderActions'
import { addInputNode } from '../FlowElement/FlowElementActions'
import {
  deleteFlowElements,
  deleteFlowElementsById,
} from '../FlowElement/FlowElementSlice'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import { isNodeData } from '../FlowElement/FlowElementUtils'
import { setInputNodeFilePath } from './InputNodeActions'
import {
  CsvInputParamType,
  FILE_TYPE_SET,
  InputNode,
  INPUT_NODE_SLICE_NAME,
} from './InputNodeType'
import { isCsvInputNode, isHDF5InputNode } from './InputNodeUtils'

const initialState: InputNode = {
  [INITIAL_IMAGE_ELEMENT_ID]: {
    fileType: FILE_TYPE_SET.IMAGE,
    param: {},
  },
}

export const inputNodeSlice = createSlice({
  name: INPUT_NODE_SLICE_NAME,
  initialState,
  reducers: {
    deleteInputNode(state, action: PayloadAction<string>) {
      delete state[action.payload]
    },
    setCsvInputNodeParam(
      state,
      action: PayloadAction<{
        nodeId: string
        param: CsvInputParamType
      }>,
    ) {
      const { nodeId, param } = action.payload
      const inputNode = state[nodeId]
      if (isCsvInputNode(inputNode)) {
        inputNode.param = param
      }
    },
    setInputNodeHDF5Path(
      state,
      action: PayloadAction<{
        nodeId: string
        path: string
      }>,
    ) {
      const { nodeId, path } = action.payload
      const item = state[nodeId]
      if (isHDF5InputNode(item)) {
        item.hdf5Path = path
      }
    },
  },
  extraReducers: (builder) =>
    builder
      .addCase(setInputNodeFilePath, (state, action) => {
        const { nodeId, filePath } = action.payload
        const targetNode = state[nodeId]
        targetNode.selectedFilePath = filePath
        if (isHDF5InputNode(targetNode)) {
          targetNode.hdf5Path = undefined
        }
      })
      .addCase(addInputNode, (state, action) => {
        const { node, fileType } = action.payload
        if (node.data?.type === NODE_TYPE_SET.INPUT) {
          switch (fileType) {
            case FILE_TYPE_SET.CSV:
              state[node.id] = {
                fileType,
                param: {
                  setHeader: null,
                  setIndex: false,
                  transpose: false,
                },
              }
              break
            case FILE_TYPE_SET.IMAGE:
              state[node.id] = {
                fileType,
                param: {},
              }
              break
            case FILE_TYPE_SET.HDF5:
              state[node.id] = {
                fileType,
                param: {},
              }
              break
            case FILE_TYPE_SET.FLUO:
              state[node.id] = {
                fileType: FILE_TYPE_SET.CSV,
                param: {
                  setHeader: null,
                  setIndex: false,
                  transpose: false,
                },
              }
              break
            case FILE_TYPE_SET.BEHAVIOR:
              state[node.id] = {
                fileType: FILE_TYPE_SET.CSV,
                param: {
                  setHeader: null,
                  setIndex: false,
                  transpose: false,
                },
              }
              break
          }
        }
      })
      .addCase(deleteFlowElements, (state, action) => {
        action.payload
          .filter((node) => isNodeData(node))
          .forEach((node) => {
            if (node.data?.type === NODE_TYPE_SET.INPUT) {
              delete state[node.id]
            }
          })
      })
      .addCase(deleteFlowElementsById, (state, action) => {
        if (Object.keys(state).includes(action.payload)) {
          delete state[action.payload]
        }
      })
      .addCase(importExperimentByUid.fulfilled, (_, action) => {
        const newState: InputNode = {}
        Object.values(action.payload.nodeDict)
          .filter(isInputNodePostData)
          .forEach((node) => {
            if (node.data != null) {
              if (node.data.fileType === FILE_TYPE_SET.IMAGE) {
                newState[node.id] = {
                  fileType: FILE_TYPE_SET.IMAGE,
                  selectedFilePath: node.data.path as string[],
                  param: {},
                }
              } else if (node.data.fileType === FILE_TYPE_SET.CSV) {
                newState[node.id] = {
                  fileType: FILE_TYPE_SET.CSV,
                  selectedFilePath: node.data.path as string,
                  param: node.data.param as CsvInputParamType,
                }
              } else if (node.data.fileType === FILE_TYPE_SET.HDF5) {
                newState[node.id] = {
                  fileType: FILE_TYPE_SET.HDF5,
                  hdf5Path: node.data.hdf5Path,
                  selectedFilePath: node.data.path as string,
                  param: {},
                }
              }
            }
          })
        return newState
      })
      .addCase(fetchExperiment.fulfilled, (_, action) => {
        const newState: InputNode = {}
        Object.values(action.payload.nodeDict)
          .filter(isInputNodePostData)
          .forEach((node) => {
            if (node.data != null) {
              if (node.data.fileType === FILE_TYPE_SET.IMAGE) {
                newState[node.id] = {
                  fileType: FILE_TYPE_SET.IMAGE,
                  selectedFilePath: node.data.path as string[],
                  param: {},
                }
              } else if (node.data.fileType === FILE_TYPE_SET.CSV) {
                newState[node.id] = {
                  fileType: FILE_TYPE_SET.CSV,
                  selectedFilePath: node.data.path as string,
                  param: node.data.param as CsvInputParamType,
                }
              } else if (node.data.fileType === FILE_TYPE_SET.HDF5) {
                newState[node.id] = {
                  fileType: FILE_TYPE_SET.HDF5,
                  hdf5Path: node.data.hdf5Path,
                  selectedFilePath: node.data.path as string,
                  param: {},
                }
              }
            }
          })
        return newState
      })
      .addCase(uploadFile.fulfilled, (state, action) => {
        const { nodeId } = action.meta.arg
        if (nodeId != null) {
          const { resultPath } = action.payload
          const target = state[nodeId]
          if (target.fileType === FILE_TYPE_SET.IMAGE) {
            target.selectedFilePath = [resultPath]
          } else {
            target.selectedFilePath = resultPath
          }
        }
      }),
})

export const { setCsvInputNodeParam, setInputNodeHDF5Path } =
  inputNodeSlice.actions

export default inputNodeSlice.reducer
