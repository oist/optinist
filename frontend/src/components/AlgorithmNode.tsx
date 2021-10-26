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
  imagePathMaxIndexByIdSelector,
  outputKeyListSelector,
  selectedOutputKeySelector,
  selectedOutputPathTypeSelector,
  selectedOutputPathValueSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { setSelectedOutputKey } from 'redux/slice/Algorithm/Algorithm'
import { showAlgoOutputImage } from 'redux/slice/ImageIndex/ImageIndex'
import { FlexLayoutModelContext } from 'App'

import { OUTPUT_TABSET_ID, PARAM_FORM_TABSET_ID } from 'const/flexlayout'
import { useTabAction } from 'FlexLayoutHook'

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

export const AlgorithmNode = React.memo<NodeProps<NodeData>>((element) => {
  const nodeId = element.id
  const theme = useTheme()
  const dispatch = useDispatch()
  const model = React.useContext(FlexLayoutModelContext)

  const selectedOutputType = useSelector(selectedOutputPathTypeSelector(nodeId))
  const selectedOutputKey = useSelector(selectedOutputKeySelector(nodeId))
  // todo imagesとfluoで決め打ちで無くなったら改修する
  const actionForOutputTab = useTabAction(
    nodeId,
    selectedOutputType === 'image' ? 'image' : 'output',
    OUTPUT_TABSET_ID,
    selectedOutputKey ?? '',
  )

  const actionForParamFormTab = useTabAction(
    nodeId,
    'paramForm',
    PARAM_FORM_TABSET_ID,
  )
  const selectedPathValue = useSelector(selectedOutputPathValueSelector(nodeId))
  const maxIndex = useSelector(imagePathMaxIndexByIdSelector(nodeId, 'images'))
  const onClick = () => {
    if (
      selectedOutputKey != null &&
      selectedPathValue != null &&
      selectedOutputType === 'image' &&
      maxIndex != null
    ) {
      dispatch(
        showAlgoOutputImage({
          id: nodeId,
          folder: selectedPathValue,
          algoName: element.data.label,
          maxIndex,
        }),
      )
    }
    if (actionForParamFormTab != null) {
      model.doAction(actionForParamFormTab)
    }
    // selectedOutputKeyがnullの場合は実行結果が存在しないためoutput用のタブを選択or作成しない
    if (
      selectedOutputKey != null &&
      selectedOutputType != null &&
      actionForOutputTab != null
    ) {
      model.doAction(actionForOutputTab)
    }
  }

  React.useEffect(() => {
    if (
      selectedOutputKey != null &&
      selectedOutputType != null &&
      actionForOutputTab != null
    ) {
      model.doAction(actionForOutputTab)
    }
  }, [selectedOutputKey])

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
        id={nodeId + '-left'}
        style={leftHandleStyle}
      />
      <Handle
        type="source"
        position={Position.Right}
        id={nodeId + '-right'}
        style={rightHandleStyle}
      />
      <OutputKeySelect nodeId={nodeId} />
    </div>
  )
})

const OutputKeySelect = React.memo<{ nodeId: string }>(({ nodeId }) => {
  const dispatch = useDispatch()
  const outputKeyList = useSelector(outputKeyListSelector(nodeId))
  const selectedOutputKey = useSelector(selectedOutputKeySelector(nodeId))
  const handleChange = (event: React.ChangeEvent<{ value: unknown }>) => {
    const outputKey = event.target.value as string
    dispatch(
      setSelectedOutputKey({
        id: nodeId,
        outputKey,
      }),
    )
  }
  if (outputKeyList.length !== 0 && selectedOutputKey != null) {
    return (
      <Select value={selectedOutputKey} label="output" onChange={handleChange}>
        {outputKeyList.map((outputKey, i) => (
          <MenuItem value={outputKey} key={i.toFixed()}>
            {outputKey}
          </MenuItem>
        ))}
      </Select>
    )
  } else {
    return null
  }
})
