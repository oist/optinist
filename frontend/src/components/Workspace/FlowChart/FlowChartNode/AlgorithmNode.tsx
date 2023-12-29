import { memo, useContext, useRef, useState } from "react"
import { useDispatch, useSelector } from "react-redux"
import { Handle, Position, NodeProps } from "reactflow"

import CheckCircleRoundedIcon from "@mui/icons-material/CheckCircleRounded"
import ErrorIcon from "@mui/icons-material/Error"
import {
  Typography,
  useTheme,
  Tooltip,
  IconButton,
  Button,
  LinearProgress,
  ButtonGroup,
  Grid,
} from "@mui/material"

import { AlgorithmInfo } from "api/algolist/AlgoList"
import { DialogContext } from "components/Workspace/FlowChart/Dialog/DialogContext"
import {
  toHandleId,
  isValidConnection,
} from "components/Workspace/FlowChart/FlowChartNode/FlowChartUtils"
import { useHandleColor } from "components/Workspace/FlowChart/FlowChartNode/HandleColorHook"
import { NodeContainer } from "components/Workspace/FlowChart/FlowChartNode/NodeContainer"
import { HANDLE_STYLE } from "const/flowchart"
import {
  selectAlgoArgs,
  selectAlgoReturns,
} from "store/slice/AlgorithmList/AlgorithmListSelectors"
import {
  selectAlgorithmIsUpdated,
  selectAlgorithmNodeDefined,
} from "store/slice/AlgorithmNode/AlgorithmNodeSelectors"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { NodeData, NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import {
  selectPipelineLatestUid,
  selectPipelineNodeResultMessage,
  selectPipelineNodeResultStatus,
  selectPipelineStatus,
} from "store/slice/Pipeline/PipelineSelectors"
import {
  NODE_RESULT_STATUS,
  RUN_STATUS,
} from "store/slice/Pipeline/PipelineType"
import { toggleParamForm } from "store/slice/RightDrawer/RightDrawerSlice"
import { RootState } from "store/store"

export const AlgorithmNode = memo(function AlgorithmNode(
  element: NodeProps<NodeData>,
) {
  const defined = useSelector(selectAlgorithmNodeDefined(element.id))
  if (defined) {
    return <AlgorithmNodeImple {...element} />
  } else {
    return null
  }
})

const AlgorithmNodeImple = memo(function AlgorithmNodeImple({
  id: nodeId,
  selected: elementSelected,
  isConnectable,
  data,
}: NodeProps<NodeData>) {
  const { onOpenOutputDialog } = useContext(DialogContext)
  const dispatch = useDispatch()

  const onClickParamButton = () => {
    dispatch(toggleParamForm(nodeId))
  }

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowNodeById(nodeId))
  }

  const onClickOutputButton = () => {
    onOpenOutputDialog(nodeId)
  }

  const status = useStatus(nodeId)
  const workflowId = useSelector(selectPipelineLatestUid)
  const isUpdated = useSelector(selectAlgorithmIsUpdated(nodeId))
  const updated = typeof workflowId !== "undefined" && isUpdated

  return (
    <NodeContainer
      nodeId={nodeId}
      selected={elementSelected}
      updated={!!updated}
    >
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: "black", position: "absolute", top: -10, right: 10 }}
      >
        ×
      </button>
      <AlgoProgress nodeId={nodeId} />
      <Grid container paddingBottom={1} justifyContent="space-between">
        <Grid item xs={10}>
          <AlgoName nodeId={nodeId} data={data} />
        </Grid>
        <Grid item xs={2}>
          <Message nodeId={nodeId} />
        </Grid>
      </Grid>
      <ButtonGroup>
        <Button
          size="small"
          onClick={onClickParamButton}
          disabled={status === NODE_RESULT_STATUS.PENDING}
        >
          Param
        </Button>
        <Button
          size="small"
          onClick={onClickOutputButton}
          disabled={status === NODE_RESULT_STATUS.PENDING}
        >
          Output
        </Button>
      </ButtonGroup>
      <AlgoArgs nodeId={nodeId} />
      <AlgoReturns nodeId={nodeId} isConnectable={isConnectable} />
    </NodeContainer>
  )
})

const AlgoProgress = memo(function AlgoProgress({ nodeId }: NodeIdProps) {
  const status = useStatus(nodeId)
  const pipelineStatus = useSelector(selectPipelineStatus)
  if (
    pipelineStatus === RUN_STATUS.START_SUCCESS &&
    status === NODE_RESULT_STATUS.PENDING
  ) {
    return (
      <div style={{ paddingLeft: 8, paddingRight: 8 }}>
        <LinearProgress />
      </div>
    )
  } else {
    return null
  }
})

interface AlgoNameProps extends NodeIdProps {
  data: NodeData
}

const AlgoName = function AlgoName({ nodeId, data }: AlgoNameProps) {
  const theme = useTheme()
  const status = useStatus(nodeId)
  return (
    <div className="algoName">
      <Typography
        style={{
          textAlign: "left",
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
}

const AlgoArgs = memo(function AlgoArgs({ nodeId }: NodeIdProps) {
  const algoArgs = useSelector(selectAlgoArgs(nodeId), algoInfoListEqualtyFn)

  return (
    <>
      {algoArgs != null
        ? algoArgs
            .filter((info) => info.type !== "params")
            .map((algoInfo, i) => {
              return (
                <ArgHandle key={i} algoInfo={algoInfo} i={i} nodeId={nodeId} />
              )
            })
        : null}
    </>
  )
})

interface AlgoReturnsProps extends NodeIdProps {
  isConnectable: boolean
}

const AlgoReturns = memo(function AlgoReturns({
  nodeId,
  isConnectable,
}: AlgoReturnsProps) {
  const algoReturns = useSelector(
    selectAlgoReturns(nodeId),
    algoInfoListEqualtyFn,
  )
  return (
    <>
      {algoReturns != null ? (
        algoReturns?.map((algoInfo, i) => {
          return (
            <ReturnHandle key={i} algoInfo={algoInfo} i={i} nodeId={nodeId} />
          )
        })
      ) : (
        // algoReturns.lengthが0の場合の応急処置
        <Handle
          type="source"
          position={Position.Right}
          id={`${nodeId}`}
          style={{
            ...HANDLE_STYLE,
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
    const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex)
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

const ArgHandle = memo(function ArgHandle({
  algoInfo: { name, type, isNone },
  nodeId,
  i,
}: HandleProps) {
  const hex_color = useHandleColor(type)
  const id = toHandleId(nodeId, name, type)
  const [isHover, setHover] = useState(false)
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
        ...HANDLE_STYLE,
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
})

const ReturnHandle = memo<HandleProps>(function ReturnHandle({
  algoInfo: { name, type },
  nodeId,
  i,
}: HandleProps) {
  const color = useHandleColor(type)
  const id = toHandleId(nodeId, name, type)
  const [isHover, setHover] = useState(false)
  return (
    <Handle
      onMouseEnter={() => setHover(true)}
      onMouseLeave={() => setHover(false)}
      key={i.toFixed()}
      type="source"
      position={Position.Right}
      id={id}
      style={{
        ...HANDLE_STYLE,
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
})

const Message = memo(function Message({ nodeId }: NodeIdProps) {
  const status = useStatus(nodeId)
  const latestUid = useSelector(selectPipelineLatestUid)
  const errorMsg = useSelector((state: RootState) =>
    latestUid != null
      ? selectPipelineNodeResultMessage(nodeId)(state) ?? null
      : null,
  )

  const anchorElRef = useRef<HTMLButtonElement | null>(null)
  const theme = useTheme()
  const { onMessageError } = useContext(DialogContext)

  if (status === NODE_RESULT_STATUS.ERROR) {
    return (
      <IconButton
        ref={anchorElRef}
        onClick={() => {
          onMessageError({ anchorElRef, message: errorMsg as string })
        }}
        size="small"
        style={{ color: theme.palette.error.main, padding: 0 }}
      >
        <ErrorIcon />
      </IconButton>
    )
  } else if (status === NODE_RESULT_STATUS.SUCCESS) {
    return (
      <CheckCircleRoundedIcon
        color="success"
        style={{ verticalAlign: "middle" }}
      />
    )
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
      : "uninitialized",
  )
  return status
}
