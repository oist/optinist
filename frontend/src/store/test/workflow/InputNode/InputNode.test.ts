import { store, rootReducer } from 'store/store'
import { REACT_FLOW_NODE_TYPE_KEY } from 'const/flowchart'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import { addInputNode } from 'store/slice/FlowElement/FlowElementActions'
import { NODE_TYPE_SET } from 'store/slice/FlowElement/FlowElementType'
import { selectNodeById } from 'store/slice/FlowElement/FlowElementSelectors'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { selectInputNodeById } from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'

describe('InputNode', () => {
  const initialRootState = store.getState()

  const nodeId = 'input_cxg72pynkd'
  const nodeType = REACT_FLOW_NODE_TYPE_KEY.ImageFileNode
  const label = 'image node'
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
  test(deleteFlowElementsById.type, () => {
    const deleteFlowElementsByIdAction = deleteFlowElementsById(nodeId)
    const targetState = rootReducer(
      rootReducer(initialRootState, addImageInputNodeAction),
      deleteFlowElementsByIdAction,
    )
    // FlowElement
    expect(selectNodeById(nodeId)(targetState)).toBeUndefined()
    // InputNode
    expect(selectInputNodeById(nodeId)(targetState)).toBeUndefined()
  })

  // filePath変更
  test(`${setInputNodeFilePath.type} by image`, () => {
    const filePath = [
      '/tmp/optinist/input/hoge/hoge.tif',
      '/tmp/optinist/input/copy_image1/copy_image1.tif',
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
