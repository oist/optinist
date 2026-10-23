import { memo } from "react"
import { useSelector, useDispatch } from "react-redux"
import { Handle, Position, NodeProps } from "reactflow"

import { FileSelect } from "components/Workspace/FlowChart/FlowChartNode/FileSelect"
import {
  toHandleId,
  isValidConnection,
} from "components/Workspace/FlowChart/FlowChartNode/FlowChartUtils"
import { useHandleColor } from "components/Workspace/FlowChart/FlowChartNode/HandleColorHook"
import { NodeContainer } from "components/Workspace/FlowChart/FlowChartNode/NodeContainer"
import { HANDLE_STYLE } from "const/flowchart"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import {
  selectImageInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from "store/slice/InputNode/InputNodeSelectors"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { arrayEqualityFn } from "utils/EqualityUtils"

export const ImageFileNode = memo(function ImageFileNode(element: NodeProps) {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <ImageFileNodeImple {...element} />
  } else {
    return null
  }
})

const ImageFileNodeImple = memo(function ImageFileNodeImple({
  id: nodeId,
  selected: elementSelected,
}: NodeProps) {
  const dispatch = useDispatch()
  const filePath = useSelector(
    selectImageInputNodeSelectedFilePath(nodeId),
    (a, b) =>
      a != null && b != null
        ? arrayEqualityFn(a as string[], b as string[])
        : a === b,
  )
  const onChangeFilePath = (path: string[]) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }

  const returnType = "ImageData"
  const imageColor = useHandleColor(returnType)

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowNodeById(nodeId))
  }

  return (
    <NodeContainer nodeId={nodeId} selected={elementSelected}>
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: "black", position: "absolute", top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        nodeId={nodeId}
        multiSelect
        onChangeFilePath={(path) => {
          if (Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        fileType={FILE_TYPE_SET.IMAGE}
        filePath={typeof filePath === "string" ? [filePath] : filePath || []}
      />
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, "image", returnType)}
        style={{
          ...HANDLE_STYLE,
          background: imageColor,
        }}
        isValidConnection={isValidConnection}
      />
    </NodeContainer>
  )
})
