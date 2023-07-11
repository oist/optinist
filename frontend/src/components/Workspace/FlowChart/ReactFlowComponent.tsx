import React, { DragEvent, MouseEvent, useState } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import ReactFlow, {
  ReactFlowProvider,
  addEdge,
  Controls,
  Connection,
  Edge,
  Node,
  OnNodesChange,
  OnEdgesChange,
  ReactFlowInstance,
  OnMove,
  Viewport,
} from 'react-flow-renderer'
import { useDrop } from 'react-dnd'

import 'style/flow.css'
import {
  editFlowNodePositionById,
  setEdgesChange,
  setFlowEdges,
  setFlowPosition,
  setNodesChange,
} from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectFlowEdges,
  selectFlowNodes,
  selectFlowPosition,
} from 'store/slice/FlowElement/FlowElementSelectors'
import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { ToolBar } from 'components/Workspace/ToolBar'
import {
  reactFlowEdgeTypes,
  reactFlowNodeTypes,
} from './FlowChartNode/ReactFlowNodeTypesConst'
import {
  DND_ITEM_TYPE_SET,
  TreeItemCollectedProps,
  TreeItemDragObject,
  TreeItemDropResult,
} from './DnDItemType'
import { AlgorithmOutputDialog } from './FlowChartNode/AlgorithmOutputDialog'
import {
  DialogContext,
  ErrorDialogValue,
  OpenDialogValue,
} from 'components/Workspace/FlowChart/DialogContext'
import { FileSelectDialog } from 'components/common/FileSelectDialog'
import { FormHelperText, Popover } from '@mui/material'

const initDialogFile = {
  filePath: '',
  open: false,
  fileTreeType: undefined,
  multiSelect: false,
  onSelectFile: () => null,
}

export const ReactFlowComponent = React.memo<UseRunPipelineReturnType>(
  (props) => {
    const flowNodes = useSelector(selectFlowNodes)
    const flowEdges = useSelector(selectFlowEdges)
    const dispatch = useDispatch()
    const [dialogNodeId, setDialogNodeId] = useState('')
    const [dialogFile, setDialogFile] =
      useState<OpenDialogValue>(initDialogFile)
    const [messageError, setMessageError] = useState<ErrorDialogValue>({
      anchorElRef: { current: null },
      message: '',
    })

    const onConnect = (params: Connection | Edge) => {
      dispatch(
        setFlowEdges(
          addEdge(
            {
              ...params,
              animated: false,
              style: { width: 5 },
              type: 'buttonedge',
            },
            flowEdges,
          ),
        ),
      )
    }

    const onNodesChange: OnNodesChange = (changes) =>
      dispatch(setNodesChange(changes))
    const onEdgesChange: OnEdgesChange = (changes) =>
      dispatch(setEdgesChange(changes))

    const onDragOver = (event: DragEvent) => {
      event.preventDefault()
      event.dataTransfer.dropEffect = 'move'
    }

    const onNodeDragStop = (event: MouseEvent, node: Node) => {
      dispatch(
        editFlowNodePositionById({
          nodeId: node.id,
          coord: { x: node.position.x, y: node.position.y },
        }),
      )
    }

    const flowPosition = useSelector(selectFlowPosition)

    const onMoveEnd: OnMove = (event, viewport: Viewport) => {
      if (event !== undefined) {
        dispatch(setFlowPosition([viewport.x, viewport.y, viewport.zoom]))
      }
    }

    const [reactFlowInstance, setReactFlowInstance] =
      React.useState<ReactFlowInstance>()

    const onInit = (reactFlowInstance: ReactFlowInstance) =>
      setReactFlowInstance(reactFlowInstance)

    const wrapparRef = React.useRef<HTMLDivElement>(null)
    const [, drop] = useDrop<
      TreeItemDragObject,
      TreeItemDropResult,
      TreeItemCollectedProps
    >(
      () => ({
        accept: DND_ITEM_TYPE_SET.TREE_ITEM,
        drop: (_, monitor) => {
          let position: TreeItemDropResult['position'] = undefined
          const monitorOffset = monitor.getClientOffset()
          if (
            wrapparRef.current != null &&
            monitorOffset != null &&
            reactFlowInstance != null
          ) {
            position = reactFlowInstance.project({
              x: monitorOffset.x - wrapparRef.current.offsetLeft - 40,
              y: monitorOffset.y - wrapparRef.current.offsetTop - 40,
            })
          }
          return { position }
        },
      }),
      [reactFlowInstance],
    )
    return (
      <div className="flow">
        <DialogContext.Provider
          value={{
            onOpen: setDialogNodeId,
            onOpenDialogFile: setDialogFile,
            onMessageError: setMessageError,
          }}
        >
          <ReactFlowProvider>
            <div className="reactflow-wrapper" ref={wrapparRef}>
              <ReactFlow
                ref={drop}
                nodes={flowNodes}
                edges={flowEdges}
                onNodesChange={onNodesChange}
                onEdgesChange={onEdgesChange}
                onConnect={onConnect}
                onInit={onInit}
                onDragOver={onDragOver}
                onNodeDragStop={onNodeDragStop}
                nodeTypes={reactFlowNodeTypes}
                edgeTypes={reactFlowEdgeTypes}
                defaultPosition={[flowPosition[0], flowPosition[1]]}
                defaultZoom={flowPosition[2]}
                onMoveEnd={onMoveEnd}
              >
                <ToolBar {...props} />
                <Controls />
              </ReactFlow>
            </div>
          </ReactFlowProvider>
          {dialogNodeId && (
            <AlgorithmOutputDialog
              nodeId={dialogNodeId}
              open
              onClose={() => setDialogNodeId('')}
            />
          )}
          {dialogFile.open && (
            <FileSelectDialog
              multiSelect={dialogFile.multiSelect}
              initialFilePath={dialogFile.filePath}
              open={dialogFile.open}
              onClickOk={(path) => {
                dialogFile.onSelectFile(path)
                setDialogFile(initDialogFile)
              }}
              onClickCancel={() => {
                setDialogFile(initDialogFile)
              }}
              fileType={dialogFile.fileTreeType}
            />
          )}
          {messageError?.message && (
            <Popover
              open
              anchorEl={messageError.anchorElRef.current}
              onClose={() =>
                setMessageError({
                  anchorElRef: { current: null },
                  message: '',
                })
              }
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
                <FormHelperText error={true}>
                  {messageError.message}
                </FormHelperText>
              </div>
            </Popover>
          )}
        </DialogContext.Provider>
      </div>
    )
  },
)
