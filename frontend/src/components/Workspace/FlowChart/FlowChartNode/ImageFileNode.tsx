import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { Handle, Position, NodeProps } from 'reactflow'

import {
  selectImageInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from 'store/slice/InputNode/InputNodeSelectors'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'

import { useHandleColor } from './HandleColorHook'
import { FileSelect } from './FileSelect'
import { toHandleId, isValidConnection } from './FlowChartUtils'
import { deleteFlowNodeById } from 'store/slice/FlowElement/FlowElementSlice'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { HANDLE_STYLE } from 'const/flowchart'
import { NodeContainer } from 'components/Workspace/FlowChart/FlowChartNode/NodeContainer'

export const ImageFileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <ImageFileNodeImple {...element} />
  } else {
    return null
  }
})

const ImageFileNodeImple = React.memo<NodeProps>(
  ({ id: nodeId, selected: elementSelected }) => {
    const dispatch = useDispatch()
    const filePath = useSelector(
      selectImageInputNodeSelectedFilePath(nodeId),
      (a, b) => (a != null && b != null ? arrayEqualityFn(a, b) : a === b),
    )
    const onChangeFilePath = (path: string[]) => {
      dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
    }

    const returnType = 'ImageData'
    const imageColor = useHandleColor(returnType)

    const onClickDeleteIcon = () => {
      dispatch(deleteFlowNodeById(nodeId))
    }

    return (
      <NodeContainer nodeId={nodeId} selected={elementSelected}>
        <button
          className="flowbutton"
          onClick={onClickDeleteIcon}
          style={{ color: 'black', position: 'absolute', top: -10, right: 10 }}
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
          filePath={filePath ?? []}
        />
        <Handle
          type="source"
          position={Position.Right}
          id={toHandleId(nodeId, 'image', returnType)}
          style={{
            ...HANDLE_STYLE,
            background: imageColor,
          }}
          isValidConnection={isValidConnection}
        />
      </NodeContainer>
    )
  },
)
