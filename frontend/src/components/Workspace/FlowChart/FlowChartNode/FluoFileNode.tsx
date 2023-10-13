import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'reactflow'

import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectCsvInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import { deleteFlowNodeById } from 'store/slice/FlowElement/FlowElementSlice'
import { useHandleColor } from './HandleColorHook'
import { ParamSettingDialog } from './CsvFileNode'
import { HANDLE_STYLE } from 'const/flowchart'
import { NodeContainer } from 'components/Workspace/FlowChart/FlowChartNode/NodeContainer'

export const FluoFileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <FluoFileNodeImple {...element} />
  } else {
    return null
  }
})

const FluoFileNodeImple = React.memo<NodeProps>(({ id: nodeId, selected }) => {
  const dispatch = useDispatch()
  const filePath = useSelector(selectCsvInputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }
  const returnType = 'FluoData'
  const fluoColor = useHandleColor(returnType)

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowNodeById(nodeId))
  }

  return (
    <NodeContainer selected={selected}>
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: 'black', position: 'absolute', top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        nodeId={nodeId}
        onChangeFilePath={(path) => {
          if (!Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        fileType={FILE_TYPE_SET.CSV}
        filePath={filePath ?? ''}
      />
      {!!filePath && <ParamSettingDialog nodeId={nodeId} filePath={filePath} />}
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, 'fluo', returnType)}
        style={{
          ...HANDLE_STYLE,
          background: fluoColor,
        }}
      />
    </NodeContainer>
  )
})
