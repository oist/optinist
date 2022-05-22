import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import {
  alpha,
  Typography,
  useTheme,
  Tooltip,
  FormHelperText,
  IconButton,
  Button,
  LinearProgress,
  ButtonGroup,
} from '@mui/material'
import CheckCircleRoundedIcon from '@mui/icons-material/CheckCircleRounded'
import ErrorIcon from '@mui/icons-material/Error'
import Popover from '@mui/material/Popover'

import { AlgorithmInfo } from 'api/algolist/AlgoList'
import {
  selectAlgoArgs,
  selectAlgoReturns,
} from 'store/slice/AlgorithmList/AlgorithmListSelectors'
import { selectAlgorithmNodeDefined } from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { NodeData } from 'store/slice/FlowElement/FlowElementType'

import { useHandleColor } from './HandleColorHook'
import { toHandleId, isValidConnection } from './FlowChartUtils'
import { toggleParamForm } from 'store/slice/RightDrawer/RightDrawerSlice'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectPipelineLatestUid,
  selectPipelineNodeResultMessage,
  selectPipelineNodeResultStatus,
  selectPipelineStatus,
} from 'store/slice/Pipeline/PipelineSelectors'
import { RootState } from 'store/store'
import {
  NODE_RESULT_STATUS,
  RUN_STATUS,
} from 'store/slice/Pipeline/PipelineType'
import { AlgorithmOutputDialog } from './AlgorithmOutputDialog'

const leftHandleStyle: CSSProperties = {
  width: '4%',
  height: '13%',
  border: '1px solid',
  borderRadius: 0,
}
const rightHandleStyle: CSSProperties = {
  width: '4%',
  height: '13%',
  border: '1px solid',
  borderRadius: 0,
}

export const AlgorithmNode = React.memo<NodeProps<NodeData>>((element) => {
  const defined = useSelector(selectAlgorithmNodeDefined(element.id))
  if (defined) {
    return <AlgorithmNodeImple {...element} />
  } else {
    return null
  }
})

const AlgorithmNodeImple = React.memo<NodeProps<NodeData>>(
  ({ id: nodeId, selected: elementSelected, isConnectable, data }) => {
    const theme = useTheme()
    const dispatch = useDispatch()

    const onClickParamButton = () => {
      dispatch(toggleParamForm(nodeId))
    }

    const onClickDeleteIcon = () => {
      dispatch(deleteFlowElementsById(nodeId))
    }

    const [open, setOpen] = React.useState(false)
    const onCloseOutputDialog = () => {
      setOpen(false)
    }
    const onClickOutputButton = () => {
      setOpen(true)
    }

    const status = useStatus(nodeId)

    return (
      <div
        style={{
          width: '100%',
          height: '110%',
          background: elementSelected
            ? alpha(theme.palette.primary.light, 0.15)
            : undefined,
          border: '1px solid',
        }}
      >
        <button
          className="flowbutton"
          onClick={onClickDeleteIcon}
          style={{ color: 'black', position: 'absolute', top: -10, right: 10 }}
        >
          ×
        </button>
        <AlgoName nodeId={nodeId} data={data} />
        <ButtonGroup sx={{ mx: 1 }}>
          <Button size="small" onClick={onClickParamButton}>
            Param
          </Button>
          <Button
            size="small"
            onClick={onClickOutputButton}
            disabled={status !== NODE_RESULT_STATUS.SUCCESS}
          >
            Output
          </Button>
          <AlgorithmOutputDialog
            nodeId={nodeId}
            open={open}
            onClose={onCloseOutputDialog}
          />
        </ButtonGroup>
        <AlgoArgs nodeId={nodeId} />
        <AlgoReturns nodeId={nodeId} isConnectable={isConnectable} />
        <Message nodeId={nodeId} />
      </div>
    )
  },
)

const AlgoName = React.memo<{
  nodeId: string
  data: NodeData
}>(({ nodeId, data }) => {
  const theme = useTheme()
  const status = useStatus(nodeId)
  const pipelineStatus = useSelector(selectPipelineStatus)
  return (
    <div
      style={{
        padding: 8,
        paddingLeft: 8,
      }}
      className="algoName"
    >
      {pipelineStatus === RUN_STATUS.START_SUCCESS &&
        status === NODE_RESULT_STATUS.PENDING && <LinearProgress />}
      <Typography
        style={{
          textAlign: 'left',
          color:
            status === NODE_RESULT_STATUS.ERROR
              ? theme.palette.error.main
              : undefined,
        }}
      >
        {data.label}
      </Typography>
    </div>
  )
})

const AlgoArgs = React.memo<{
  nodeId: string
}>(({ nodeId }) => {
  const algoArgs = useSelector(selectAlgoArgs(nodeId), algoInfoListEqualtyFn)

  return (
    <>
      {algoArgs != null
        ? algoArgs
            .filter((info) => info.type !== 'params')
            .map((algoInfo, i) => {
              return <ArgHandle algoInfo={algoInfo} i={i} nodeId={nodeId} />
            })
        : null}
    </>
  )
})

const AlgoReturns = React.memo<{
  nodeId: string
  isConnectable: boolean
}>(({ nodeId, isConnectable }) => {
  const algoReturns = useSelector(
    selectAlgoReturns(nodeId),
    algoInfoListEqualtyFn,
  )
  return (
    <>
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
    </>
  )
})

type HandleProps = {
  algoInfo: AlgorithmInfo
  nodeId: string
  i: number
}

function hexToRgb(hex: string | undefined, isNone: boolean | undefined) {
  if (hex !== undefined) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex)
    if (result !== null) {
      if (isNone) {
        return `rgba(${parseInt(result[1], 16)}, ${parseInt(
          result[2],
          16,
        )}, ${parseInt(result[3], 16)}, 0.55)`
      } else {
        return `rgba(${parseInt(result[1], 16)}, ${parseInt(
          result[2],
          16,
        )}, ${parseInt(result[3], 16)}, 1)`
      }
    } else {
      return undefined
    }
  } else {
    return undefined
  }
}

const ArgHandle = React.memo<HandleProps>(
  ({ algoInfo: { name, type, isNone }, nodeId, i }) => {
    const hex_color = useHandleColor(type)
    const id = toHandleId(nodeId, name, type)
    const [isHover, setHover] = React.useState(false)
    const rgb_color = hexToRgb(hex_color, isNone)
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
          background: rgb_color,
          top: i * 25 + 15,
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
          top: i * 25 + 15,
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

const Message = React.memo<{
  nodeId: string
}>(({ nodeId }) => {
  const status = useStatus(nodeId)
  const latestUid = useSelector(selectPipelineLatestUid)
  const errorMsg = useSelector((state: RootState) =>
    latestUid != null
      ? selectPipelineNodeResultMessage(nodeId)(state) ?? null
      : null,
  )

  const anchorElRef = React.useRef<HTMLButtonElement | null>(null)
  const [open, setOpen] = React.useState(false)
  const theme = useTheme()

  if (status === NODE_RESULT_STATUS.ERROR) {
    return (
      <>
        <IconButton
          ref={anchorElRef}
          onClick={() => setOpen((prevOpen) => !prevOpen)}
          size="small"
          style={{ color: theme.palette.error.main, float: 'right' }}
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
            <FormHelperText error={true}>{errorMsg}</FormHelperText>
          </div>
        </Popover>
      </>
    )
  } else if (status === NODE_RESULT_STATUS.SUCCESS) {
    return <CheckCircleRoundedIcon color="success" sx={{ float: 'right' }} />
  } else {
    return null
  }
})

function algoInfoListEqualtyFn(
  a: AlgorithmInfo[] | undefined,
  b: AlgorithmInfo[] | undefined,
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

function useStatus(nodeId: string) {
  const latestUid = useSelector(selectPipelineLatestUid)
  const status = useSelector((state: RootState) =>
    latestUid != null
      ? selectPipelineNodeResultStatus(nodeId)(state)
      : 'uninitialized',
  )
  return status
}
