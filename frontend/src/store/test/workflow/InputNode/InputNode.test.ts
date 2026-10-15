import { expect, describe, test } from "@jest/globals"

import { REACT_FLOW_NODE_TYPE_KEY } from "config/fileTypes.config"
import { addInputNode } from "store/slice/FlowElement/FlowElementActions"
import { selectNodeById } from "store/slice/FlowElement/FlowElementSelectors"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { NODE_TYPE_SET } from "store/slice/FlowElement/FlowElementType"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import {
  selectFilePathIsUndefined,
  selectInputNodeById,
} from "store/slice/InputNode/InputNodeSelectors"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import {
  isHDF5InputNode,
  isMatlabInputNode,
} from "store/slice/InputNode/InputNodeUtils"
import { importWorkflowConfig } from "store/slice/Workflow/WorkflowActions"
import { store, rootReducer } from "store/store"

describe("InputNode", () => {
  const initialRootState = store.getState()

  const nodeId = "input_cxg72pynkd"
  const nodeType = REACT_FLOW_NODE_TYPE_KEY.ImageFileNode
  const label = "image node"
  const nodeDataType = NODE_TYPE_SET.INPUT
  const fileType = FILE_TYPE_SET.IMAGE

  const addImageInputNodeAction = addInputNode({
    node: {
      id: nodeId,
      type: nodeType,
      data: { label, type: nodeDataType },
    },
    fileType,
  })

  // 追加
  test(addInputNode.type, () => {
    const targetState = rootReducer(initialRootState, addImageInputNodeAction)
    // FlowElement
    expect(selectNodeById(nodeId)(targetState)).toBeDefined()
    expect(selectNodeById(nodeId)(targetState)?.id).toBe(nodeId)
    expect(selectNodeById(nodeId)(targetState)?.type).toBe(nodeType)
    expect(selectNodeById(nodeId)(targetState)?.data).toEqual({
      label,
      type: nodeDataType,
    })
    // InputNode
    expect(selectInputNodeById(nodeId)(targetState).fileType).toBe(fileType)
  })

  // 削除
  test(deleteFlowNodeById.type, () => {
    const deleteFlowNodeByIdAction = deleteFlowNodeById(nodeId)
    const targetState = rootReducer(
      rootReducer(initialRootState, addImageInputNodeAction),
      deleteFlowNodeByIdAction,
    )
    // FlowElement
    expect(selectNodeById(nodeId)(targetState)).toBeUndefined()
    // InputNode
    expect(selectInputNodeById(nodeId)(targetState)).toBeUndefined()
  })

  // filePath変更
  test(`${setInputNodeFilePath.type} by image`, () => {
    const filePath = [
      "/tmp/optinist/input/hoge/hoge.tif",
      "/tmp/optinist/input/copy_image1/copy_image1.tif",
    ]
    const setInputNodeFilePathAction = setInputNodeFilePath({
      nodeId,
      filePath,
    })
    const targetState = rootReducer(
      rootReducer(initialRootState, addImageInputNodeAction),
      setInputNodeFilePathAction,
    )
    expect(selectInputNodeById(nodeId)(targetState).selectedFilePath).toEqual(
      filePath,
    )
  })
})

describe("InputNode importWorkflowConfig", () => {
  const initialRootState = store.getState()

  // backend shape: param/hdf5Path/matPath always present, null when unset
  const inputNodePostData = (
    id: string,
    fileType: string,
    data: Record<string, unknown> = {},
  ) => ({
    id,
    type: "FileNode",
    data: {
      label: id,
      type: NODE_TYPE_SET.INPUT,
      fileType,
      param: {},
      hdf5Path: null,
      matPath: null,
      ...data,
    },
    position: { x: 0, y: 0 },
  })

  const importAction = (nodeDict: Record<string, unknown>) => ({
    type: importWorkflowConfig.fulfilled.type,
    payload: { nodeDict, edgeDict: {} },
  })

  test("restores selectedFilePath, hdf5Path and matPath", () => {
    const state = rootReducer(
      initialRootState,
      importAction({
        image_array: inputNodePostData("image_array", FILE_TYPE_SET.IMAGE, {
          path: ["dir/a.tif"],
        }),
        image_string: inputNodePostData("image_string", FILE_TYPE_SET.IMAGE, {
          path: "b.tif",
        }),
        hdf5_node: inputNodePostData("hdf5_node", FILE_TYPE_SET.HDF5, {
          path: "x.nwb",
          hdf5Path: "acquisition/data",
        }),
        matlab_node: inputNodePostData("matlab_node", FILE_TYPE_SET.MATLAB, {
          path: "m.mat",
          matPath: "field1",
        }),
        csv_node: inputNodePostData("csv_node", FILE_TYPE_SET.CSV, {
          path: "c.csv",
        }),
      }),
    )

    expect(selectInputNodeById("image_array")(state).selectedFilePath).toEqual([
      "dir/a.tif",
    ])
    expect(selectInputNodeById("image_string")(state).selectedFilePath).toEqual(
      ["b.tif"],
    )
    const hdf5Node = selectInputNodeById("hdf5_node")(state)
    expect(hdf5Node.selectedFilePath).toBe("x.nwb")
    expect(isHDF5InputNode(hdf5Node) ? hdf5Node.hdf5Path : undefined).toBe(
      "acquisition/data",
    )
    const matlabNode = selectInputNodeById("matlab_node")(state)
    expect(matlabNode.selectedFilePath).toBe("m.mat")
    expect(isMatlabInputNode(matlabNode) ? matlabNode.matPath : undefined).toBe(
      "field1",
    )
    expect(selectInputNodeById("csv_node")(state).selectedFilePath).toBe(
      "c.csv",
    )
    expect(selectFilePathIsUndefined(state)).toBe(false)
  })

  test("null or empty paths stay unselected and keep the run gate closed", () => {
    const state = rootReducer(
      initialRootState,
      importAction({
        csv_null: inputNodePostData("csv_null", FILE_TYPE_SET.CSV, {
          path: null,
        }),
        image_empty: inputNodePostData("image_empty", FILE_TYPE_SET.IMAGE, {
          path: [""],
        }),
        hdf5_no_path: inputNodePostData("hdf5_no_path", FILE_TYPE_SET.HDF5, {
          path: "x.nwb",
          hdf5Path: "",
        }),
      }),
    )

    expect(
      selectInputNodeById("csv_null")(state).selectedFilePath,
    ).toBeUndefined()
    expect(
      selectInputNodeById("image_empty")(state).selectedFilePath,
    ).toBeUndefined()
    const hdf5Node = selectInputNodeById("hdf5_no_path")(state)
    expect(
      isHDF5InputNode(hdf5Node) ? hdf5Node.hdf5Path : "not hdf5",
    ).toBeUndefined()
    expect(selectFilePathIsUndefined(state)).toBe(true)
  })

  test("a matlab node without matPath keeps the run gate closed", () => {
    const state = rootReducer(
      initialRootState,
      importAction({
        matlab_node: inputNodePostData("matlab_node", FILE_TYPE_SET.MATLAB, {
          path: "m.mat",
        }),
      }),
    )
    expect(selectFilePathIsUndefined(state)).toBe(true)
  })
})
