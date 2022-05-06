import React, { CSSProperties } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@mui/material/styles'

import {
  selectImageInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from 'store/slice/InputNode/InputNodeSelectors'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'

import { useHandleColor } from './HandleColorHook'
import { FileSelect } from './FileSelect'
import { toHandleId, isValidConnection } from './FlowChartUtils'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'
import { arrayEqualityFn } from 'utils/EqualityUtils'

// Connection部分のレイアウト
const sourceHandleStyle: CSSProperties = {
  width: '4%',
  height: '13%',
  top: 15,
  border: '1px solid',
  // borderColor: '#555',
  borderRadius: 0,
}

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

    const theme = useTheme()
    const returnType = 'ImageData'
    const imageColor = useHandleColor(returnType)

    const onClickDeleteIcon = () => {
      dispatch(deleteFlowElementsById(nodeId))
    }

    return (
      <div
        style={{
          height: '100%',
          width: '250px',
          background: elementSelected
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
            ...sourceHandleStyle,
            background: imageColor,
          }}
          isValidConnection={isValidConnection}
        />
      </div>
    )
  },
)
