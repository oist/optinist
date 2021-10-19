import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { NodeData } from 'const/NodeData'

import {
  alpha,
  Typography,
  useTheme,
  Select,
  MenuItem,
} from '@material-ui/core'
import {
  imageDirMaxIndexByIdSelector,
  outputPathListSelector,
  selectedOutputPathSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { setSelectedOutputPath } from 'redux/slice/Algorithm/Algorithm'
import { showAlgoOutputImage } from 'redux/slice/ImageIndex/ImageIndex'

export const AlgorithmNode = React.memo<NodeProps<NodeData>>((element) => {
  const theme = useTheme()
  const dispatch = useDispatch()
  const leftHandleStyle: CSSProperties = {
    width: 8,
    left: 0,
    height: '100%',
    border: '1px solid',
    borderRadius: 0,
  }
  const rightHandleStyle: CSSProperties = {
    width: 8,
    right: 0,
    height: '100%',
    border: '1px solid',
    borderColor: 'black',
    borderRadius: 0,
  }

  const pathList = useSelector(outputPathListSelector(element.id))
  const selectedPath = useSelector(selectedOutputPathSelector(element.id))
  const handleChange = (event: React.ChangeEvent<{ value: unknown }>) => {
    dispatch(
      setSelectedOutputPath({
        id: element.id,
        path: event.target.value as string,
      }),
    )
  }

  const maxIndex = useSelector(imageDirMaxIndexByIdSelector(element.id))
  const onClick = () => {
    if (
      selectedPath != null &&
      selectedPath.value != null &&
      selectedPath.isImage
    ) {
      dispatch(
        showAlgoOutputImage({
          id: element.id,
          folder: selectedPath.value,
          algoName: element.data.label,
          maxIndex: maxIndex ?? 0,
        }),
      )
    }
  }

  return (
    <div
      style={{
        width: '100%',
        height: '100%',
        background: element.selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
      onClick={onClick}
    >
      <div
        style={{
          padding: 8,
          paddingLeft: 16,
        }}
      >
        <Typography style={{ textAlign: 'left' }}>
          {element.data.label}
          {/* ({element.id}) */}
        </Typography>
      </div>
      <Handle
        type="target"
        position={Position.Left}
        id={element.id + '-left'}
        style={leftHandleStyle}
      />
      <Handle
        type="source"
        position={Position.Right}
        id={element.id + '-right'}
        style={rightHandleStyle}
      />
      {pathList.length !== 0 && selectedPath != null && (
        <Select
          value={selectedPath.value}
          label="output"
          onChange={handleChange}
        >
          {pathList.map(([key, value]) => (
            <MenuItem value={typeof value === 'string' ? value : value?.path}>
              {key}
            </MenuItem>
          ))}
        </Select>
      )}
    </div>
  )
})
