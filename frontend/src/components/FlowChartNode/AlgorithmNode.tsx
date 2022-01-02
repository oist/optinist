import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
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

import { OutputPathsDTO } from 'api/Run/Run'
import { AlgorithmInfo } from 'store/slice/AlgorithmList/AlgorithmListType'
import {
  selectAlgoArgs,
  selectAlgoReturns,
} from 'store/slice/AlgorithmList/AlgorithmListSelectors'
import {
  selectAlgorithmFunctionPath,
  selectAlgorithmSelectedOutputKey,
  selectAlgorithmNodeDefined,
} from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { setSelectedOutputKey } from 'store/slice/AlgorithmNode/AlgorithmNodeSlice'
import { NodeData } from 'store/slice/FlowElement/FlowElementType'
import { toDataType } from 'store/slice/DisplayData/DisplayDataUtils'
import {
  selectOutputPaths,
  selectResultError,
} from 'store/slice/RunPipelineResult/RenPipelineResultSelectors'

import {
  useDisplayDataTabAciton,
  useParamFormTabAction,
} from 'components/flextlayout/FlexLayoutHook'
import { useHandleColor } from './HandleColorHook'
import { toHandleId, isValidConnection } from './FlowChartUtils'

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
    const outputPaths = useSelector(selectOutputPaths)
    const functionPath = useSelector(selectAlgorithmFunctionPath(nodeId))
    const selectedOutputKey = useSelector(
      selectAlgorithmSelectedOutputKey(nodeId),
    )
    const outputKeyList =
      outputPaths != null ? getOutputKeyList(outputPaths, functionPath) : null

    const setParamFormTab = useParamFormTabAction(nodeId)
    const { setDisplayTab } = useDisplayDataTabAciton(nodeId)
    React.useEffect(() => {
      if (outputPaths != null) {
        const keyList = getOutputKeyList(outputPaths, functionPath)
        const initialOutputKey = keyList.length > 0 ? keyList[0] : null
        if (selectedOutputKey === null && initialOutputKey != null) {
          // outputKeyの選択の初期値を設定
          dispatch(
            setSelectedOutputKey({
              nodeId,
              selectedOutputKey: initialOutputKey,
            }),
          )
        }
      }
    }, [
      dispatch,
      outputPaths,
      functionPath,
      setDisplayTab,
      nodeId,
      selectedOutputKey,
    ])

    React.useEffect(() => {
      // Nodeが置かれたタイミングでparamFormのタブを作成する
      setParamFormTab()
    }, [setParamFormTab])

    const selectetOutput =
      outputPaths != null && getOutputPathMap(outputPaths, functionPath) != null
        ? Object.entries(outputPaths[functionPath]).find(
            ([outputKey, outputInfo]) => outputKey === selectedOutputKey,
          )?.[1]
        : undefined
    const onClick = () => {
      setParamFormTab()
      if (selectetOutput !== undefined) {
        setDisplayTab(
          selectetOutput.path,
          toDataType(selectetOutput.type),
          selectedOutputKey ?? undefined,
        )
      }
    }

    const onChangeOutputKey = (outputKey: string) => {
      dispatch(
        setSelectedOutputKey({
          nodeId,
          selectedOutputKey: outputKey,
        }),
      )
      if (selectetOutput !== undefined) {
        setDisplayTab(
          selectetOutput.path,
          toDataType(selectetOutput.type),
          selectedOutputKey ?? selectetOutput.path,
        )
      }
    }

    const algoArgs = useSelector(selectAlgoArgs(nodeId), algoInfoListEqualtyFn)
    const algoReturns = useSelector(
      selectAlgoReturns(nodeId),
      algoInfoListEqualtyFn,
    )

    const error = useSelector(selectResultError(nodeId))
    return (
      <div
        style={{
          width: '100%',
          height: '100%',
          background: elementSelected
            ? alpha(theme.palette.primary.light, 0.15)
            : undefined,
          border: '1px solid',
        }}
        onClick={onClick}
      >
        <div
          style={{
            padding: 8,
            paddingLeft: 8,
          }}
          className="algoName"
        >
          <Typography
            style={{
              textAlign: 'left',
              color: error != null ? theme.palette.error.main : undefined,
            }}
          >
            {data.label}
            <ErrorMessage error={error} />
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
        {outputKeyList != null && outputKeyList.length > 0 && (
          <div className="outputkey">
            <OutputKeySelect
              selectedOutputKey={selectedOutputKey ?? outputKeyList[0]}
              outputKeyList={outputKeyList}
              onChange={onChangeOutputKey}
            />
          </div>
        )}
      </div>
    )
  },
)

const OutputKeySelect = React.memo<{
  outputKeyList: string[]
  selectedOutputKey: string
  onChange: (outputKey: string) => void
}>(({ outputKeyList, selectedOutputKey, onChange }) => {
  const handleChange = (event: React.ChangeEvent<{ value: unknown }>) => {
    const outputKey = event.target.value as string
    onChange(outputKey)
  }
  if (outputKeyList.length !== 0) {
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
  algoInfo: AlgorithmInfo
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

const ErrorMessage = React.memo<{
  error: string | null
}>(({ error }) => {
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

function getOutputPathMap(outputPaths: OutputPathsDTO, functionPath: string) {
  if (Object.keys(outputPaths).includes(functionPath)) {
    return outputPaths[functionPath]
  } else {
    return null
  }
}

function getOutputKeyList(outputPaths: OutputPathsDTO, functionPath: string) {
  const output = getOutputPathMap(outputPaths, functionPath)
  if (output != null) {
    return Object.keys(output)
  } else {
    return []
  }
}
