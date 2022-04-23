import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@mui/material/styles'

import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectCsvInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import { useHandleColor } from './HandleColorHook'
import { ParamSettingDialog } from './CsvFileNode'

const sourceHandleStyle: CSSProperties = {
  width: 8,
  height: 15,
  top: 15,
  border: '1px solid',
  borderColor: '#555',
  borderRadius: 0,
}

export const BehaviorFileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <BehaviorFileNodeImple {...element} />
  } else {
    return null
  }
})

const BehaviorFileNodeImple = React.memo<NodeProps>(
  ({ id: nodeId, selected }) => {
    const dispatch = useDispatch()
    const filePath = useSelector(selectCsvInputNodeSelectedFilePath(nodeId))
    const onChangeFilePath = (path: string) => {
      dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
    }
    const theme = useTheme()
    const returnType = 'BehaviorData'
    const behaviorColor = useHandleColor(returnType)

    const onClickDeleteIcon = () => {
      dispatch(deleteFlowElementsById(nodeId))
    }

    return (
      <div
        style={{
          height: '100%',
          width: '230px',
          background: selected
            ? alpha(theme.palette.primary.light, 0.1)
            : undefined,
        }}
      >
        <button
          className="flowbutton"
          onClick={onClickDeleteIcon}
          style={{ color: 'black', position: 'absolute', top: -10, right: 10 }}
        >
          ×
        </button>
        <FileSelect
          onChangeFilePath={(path) => {
            if (!Array.isArray(path)) {
              onChangeFilePath(path)
            }
          }}
          fileType={FILE_TYPE_SET.CSV}
          filePath={filePath ?? ''}
        />
        {!!filePath && (
          <ParamSettingDialog nodeId={nodeId} filePath={filePath} />
        )}
        <Handle
          type="source"
          position={Position.Right}
          id={toHandleId(nodeId, 'behavior', returnType)}
          style={{
            ...sourceHandleStyle,
            background: behaviorColor,
          }}
        />
      </div>
    )
  },
)
