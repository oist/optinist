import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps, Connection } from 'react-flow-renderer'
import { NodeData } from 'const/NodeData'
import {
  alpha,
  Typography,
  useTheme,
  Select,
  MenuItem,
  Tooltip,
  FormHelperText,
  IconButton,
} from '@material-ui/core'
import ErrorIcon from '@material-ui/icons/Error'
import Popover from '@material-ui/core/Popover'

import {
  algoArgsSelector,
  algoErrorSelector,
  algoReturnsSelector,
  outputKeyListSelector,
  selectedOutputKeySelector,
  selectedOutputPathTypeSelector,
} from 'store/slice/Algorithm/AlgorithmSelector'
import { setSelectedOutputKey } from 'store/slice/Algorithm/Algorithm'
import { FlexLayoutModelContext } from 'App'

import { arrayEqualityFn } from 'utils/EqualityUtils'
import { OUTPUT_TABSET_ID, PARAM_FORM_TABSET_ID } from 'const/flexlayout'
import { useTabAction } from 'FlexLayoutHook'
import { AlgoInfo, OUTPUT_TYPE_SET } from 'store/slice/Algorithm/AlgorithmType'
import { handleTypeColorSelector } from 'store/slice/HandleTypeColor/HandleTypeColorSelector'
import { addColor } from 'store/slice/HandleTypeColor/HandleTypeColor'
import { RootState } from 'store/store'

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
  const model = React.useContext(FlexLayoutModelContext)

  const selectedOutputType = useSelector(selectedOutputPathTypeSelector(nodeId))
  const selectedOutputKey = useSelector(selectedOutputKeySelector(nodeId))
  const actionForOutputTab = useTabAction(nodeId)

  const actionForParamFormTab = useTabAction(nodeId)

  React.useEffect(() => {
    // Nodeが置かれたタイミングでparamFormのタブを作成する
    if (actionForParamFormTab != null) {
      model.doAction(actionForParamFormTab('paramForm', PARAM_FORM_TABSET_ID))
    }
  }, [model, actionForParamFormTab])

  const onClick = () => {
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

  const algoArgs = useSelector(algoArgsSelector(nodeId), algoInfoListEqualtyFn)
  const algoReturns = useSelector(
    algoReturnsSelector(nodeId),
    algoInfoListEqualtyFn,
  )
  const isError = useSelector(
    (state: RootState) => algoErrorSelector(nodeId)(state) != null,
  )
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
        <Typography
          style={{
            textAlign: 'left',
            color: isError ? theme.palette.error.main : undefined,
          }}
        >
          {element.data.label}
          <ErrorMessage nodeId={nodeId} />
        </Typography>
      </div>
      <div>
        {algoArgs != null
          ? algoArgs
              .filter((info) => info.type !== 'params')
              .map((algoInfo, i) => {
                return <ArgHandle algoInfo={algoInfo} i={i} nodeId={nodeId} />
              })
          : null}
      </div>
      {algoReturns != null ? (
        algoReturns?.map((algoInfo, i) => {
          return <ReturnHandle algoInfo={algoInfo} i={i} nodeId={nodeId} />
        })
      ) : (
        // algoReturns.lengthが0の場合の応急処置
        <Handle
          type="source"
          position={Position.Right}
          id={`${nodeId}`}
          style={{
            ...rightHandleStyle,
            top: 15,
          }}
          isConnectable={isConnectable}
        />
      )}
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

type HandleProps = {
  algoInfo: AlgoInfo
  nodeId: string
  i: number
}

const ArgHandle = React.memo<HandleProps>(
  ({ algoInfo: { name, type }, nodeId, i }) => {
    const color = useHandleColor(type)
    const id = toHandleId(nodeId, name, type)
    const [isHover, setHover] = React.useState(false)
    return (
      <Handle
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
        key={i.toFixed()}
        type="target"
        position={Position.Left}
        id={id}
        style={{
          ...leftHandleStyle,
          background: color,
          top: i * 35 + 15,
        }}
        isValidConnection={isValidConnection}
      >
        <Tooltip
          title={
            <>
              <Typography color="inherit">name: {name}</Typography>
              <Typography color="inherit">type: {type}</Typography>
            </>
          }
          open={isHover}
          placement="left-end"
          arrow
        >
          <div />
        </Tooltip>
      </Handle>
    )
  },
)

const ReturnHandle = React.memo<HandleProps>(
  ({ algoInfo: { name, type }, nodeId, i }) => {
    const color = useHandleColor(type)
    const id = toHandleId(nodeId, name, type)
    const [isHover, setHover] = React.useState(false)
    return (
      <Handle
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
        key={i.toFixed()}
        type="source"
        position={Position.Right}
        id={id}
        style={{
          ...rightHandleStyle,
          background: color,
          top: i * 35 + 15,
        }}
        isValidConnection={isValidConnection}
      >
        <Tooltip
          title={
            <>
              <Typography color="inherit">name: {name}</Typography>
              <Typography color="inherit">type: {type}</Typography>
            </>
          }
          open={isHover}
          placement="right-end"
          arrow
        >
          <div />
        </Tooltip>
      </Handle>
    )
  },
)

const ErrorMessage = React.memo<{
  nodeId: string
}>(({ nodeId }) => {
  const error = useSelector(algoErrorSelector(nodeId))
  const anchorElRef = React.useRef<HTMLButtonElement | null>(null)
  const [open, setOpen] = React.useState(false)
  const theme = useTheme()
  if (error != null) {
    return (
      <>
        <IconButton
          ref={anchorElRef}
          onClick={() => setOpen((prevOpen) => !prevOpen)}
          size="small"
          style={{ color: theme.palette.error.main }}
        >
          <ErrorIcon />
        </IconButton>
        <Popover
          open={open}
          anchorEl={anchorElRef.current}
          onClose={() => setOpen(false)}
          anchorOrigin={{
            vertical: 'top',
            horizontal: 'right',
          }}
          transformOrigin={{
            vertical: 'bottom',
            horizontal: 'left',
          }}
        >
          <div style={{ margin: 8 }}>
            <FormHelperText error={true}>{error}</FormHelperText>
          </div>
        </Popover>
      </>
    )
  } else {
    return null
  }
})

function toHandleId(nodeId: string, name: string, type: string) {
  return `${nodeId}--${name}--${type}`
}

function getHandleType(handleId: string) {
  return handleId.split('--')[2]
}

function isValidConnection(connection: Connection) {
  if (connection.sourceHandle != null && connection.targetHandle != null) {
    return (
      getHandleType(connection.sourceHandle) ===
      getHandleType(connection.targetHandle)
    )
  } else {
    return true
  }
}

function useHandleColor(type: string) {
  const dispatch = useDispatch()
  const color = useSelector(handleTypeColorSelector(type))
  React.useEffect(() => {
    if (color === undefined) {
      dispatch(addColor(type))
    }
  }, [type, color, dispatch])
  return color
}

function algoInfoListEqualtyFn(
  a: AlgoInfo[] | undefined,
  b: AlgoInfo[] | undefined,
) {
  if (a != null && b != null) {
    return (
      a === b ||
      (a.length === b.length &&
        a.every((v, i) => v.type === b[i].type && v.name === b[i].name))
    )
  } else {
    return a === undefined && b === undefined
  }
}
