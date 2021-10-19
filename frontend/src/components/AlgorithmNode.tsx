import React, { CSSProperties } from 'react'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { NodeData } from 'const/NodeData'

import {
  alpha,
  Typography,
  useTheme,
  Select,
  MenuItem,
} from '@material-ui/core'

export const AlgorithmNode = React.memo<NodeProps<NodeData>>((element) => {
  const theme = useTheme()
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

  const [age, setAge] = React.useState(1)
  const handleChange = (event: any) => {
    setAge(event.target.value)
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
      <Select
        labelId="demo-simple-select-label"
        id="demo-simple-select"
        value={age}
        label="Age"
        onChange={handleChange}
      >
        <MenuItem value={1}>Images</MenuItem>
        <MenuItem value={2}>Path</MenuItem>
        <MenuItem value={3}>Result</MenuItem>
      </Select>
    </div>
  )
})
