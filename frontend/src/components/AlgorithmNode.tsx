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
  algoArgsSelector,
  algoReturnsSelector,
  imagePathMaxIndexByIdSelector,
  outputKeyListSelector,
  selectedOutputKeySelector,
  selectedOutputPathTypeSelector,
  selectedOutputPathValueSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { setSelectedOutputKey } from 'redux/slice/Algorithm/Algorithm'
import { showAlgoOutputImage } from 'redux/slice/ImageIndex/ImageIndex'
import { FlexLayoutModelContext } from 'App'

import { arrayEqualityFn } from 'utils/EqualityUtils'
import { OUTPUT_TABSET_ID, PARAM_FORM_TABSET_ID } from 'const/flexlayout'
import { useTabAction } from 'FlexLayoutHook'
import { OUTPUT_TYPE_SET } from 'redux/slice/Algorithm/AlgorithmType'

const leftHandleStyle: CSSProperties = {
  width: 8,
  height: '15%',
  border: '1px solid',
  borderRadius: 0,
}
const rightHandleStyle: CSSProperties = {
  width: 8,
  height: '15%',
  border: '1px solid',
  borderColor: 'black',
  borderRadius: 0,
}

export const AlgorithmNode = React.memo<NodeProps<NodeData>>((element) => {
  const nodeId = element.id
  const { isConnectable } = element
  const theme = useTheme()
  const dispatch = useDispatch()
  const model = React.useContext(FlexLayoutModelContext)

  const selectedOutputType = useSelector(selectedOutputPathTypeSelector(nodeId))
  const selectedOutputKey = useSelector(selectedOutputKeySelector(nodeId))
  const actionForOutputTab = useTabAction(nodeId)

  const actionForParamFormTab = useTabAction(nodeId)
  const selectedPathValue = useSelector(selectedOutputPathValueSelector(nodeId))
  const maxIndex = useSelector(
    imagePathMaxIndexByIdSelector(nodeId, selectedOutputKey ?? ''),
  )
  const onClick = () => {
    if (
      selectedOutputKey != null &&
      selectedPathValue != null &&
      selectedOutputType === OUTPUT_TYPE_SET.IMAGE &&
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
      model.doAction(actionForParamFormTab('paramForm', PARAM_FORM_TABSET_ID))
    }
    // selectedOutputKeyがnullの場合は実行結果が存在しないためoutput用のタブを選択or作成しない
    if (
      selectedOutputKey != null &&
      selectedOutputType != null &&
      actionForOutputTab != null
    ) {
      model.doAction(
        actionForOutputTab(
          selectedOutputType === OUTPUT_TYPE_SET.IMAGE ? 'image' : 'output',
          OUTPUT_TABSET_ID,
          selectedOutputKey,
        ),
      )
    }
  }

  React.useEffect(() => {
    if (
      selectedOutputKey != null &&
      selectedOutputType != null &&
      actionForOutputTab != null
    ) {
      model.doAction(
        actionForOutputTab(
          selectedOutputType === OUTPUT_TYPE_SET.IMAGE ? 'image' : 'output',
          OUTPUT_TABSET_ID,
          selectedOutputKey,
        ),
      )
    }
  }, [selectedOutputKey, actionForOutputTab, model, selectedOutputType])

  const algoArgs = useSelector(algoArgsSelector(nodeId), (a, b) => {
    if (a != null && b != null) {
      return arrayEqualityFn(a, b)
    } else {
      return a === undefined && b === undefined
    }
  })
  const algoReturns = useSelector(algoReturnsSelector(nodeId), (a, b) => {
    if (a != null && b != null) {
      return arrayEqualityFn(a, b)
    } else {
      return a === undefined && b === undefined
    }
  })
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
      <div>
        {algoArgs != null
          ? algoArgs
              .filter((name) => name !== 'params')
              .map((argsName, i) => {
                return (
                  <Handle
                    key={i.toFixed()}
                    type="target"
                    position={Position.Left}
                    id={`${nodeId}-${argsName}`}
                    style={{
                      ...leftHandleStyle,
                      top: i * 35 + 15,
                    }}
                    isConnectable={isConnectable}
                  />
                )
              })
          : null}
      </div>
      {algoReturns != null ? (
        algoReturns.length > 0 ? (
          algoReturns.map((returnsName, i) => {
            return (
              <Handle
                key={i.toFixed()}
                type="source"
                position={Position.Right}
                id={`${nodeId}-${returnsName}`}
                style={{
                  ...rightHandleStyle,
                  top: i * 35 + 15,
                }}
                isConnectable={isConnectable}
              />
            )
          })
        ) : (
          // algoReturns.lengthが0の場合の応急処置
          <Handle
            type="source"
            position={Position.Right}
            id={`${nodeId}`}
            style={{
              ...rightHandleStyle,
              height: '100%',
            }}
            isConnectable={isConnectable}
          />
        )
      ) : null}
      <OutputKeySelect nodeId={nodeId} />
    </div>
  )
})

const OutputKeySelect = React.memo<{ nodeId: string }>(({ nodeId }) => {
  const dispatch = useDispatch()
  const outputKeyList = useSelector(
    outputKeyListSelector(nodeId),
    arrayEqualityFn,
  )
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
