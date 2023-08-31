import React, {
  DragEvent,
  MouseEvent as ReactMouseEvent,
  FC,
  ReactNode,
  useState,
} from 'react'
import { useSelector, useDispatch } from 'react-redux'
import ReactFlow, {
  ReactFlowProvider,
  Controls,
  Connection,
  Edge,
  Node,
  applyNodeChanges,
  NodeTypes,
  EdgeTypes,
  NodeChange,
  Viewport,
  isNode,
  addEdge,
  ReactFlowInstance,
} from 'reactflow'
import { useDrop } from 'react-dnd'
import 'reactflow/dist/style.css'
import 'style/flow.css'
import {
  editFlowElementPositionById,
  setFlowElements,
  setFlowPosition,
} from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectFlowElements,
  selectFlowPosition,
} from 'store/slice/FlowElement/FlowElementSelectors'
import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { ToolBar } from 'components/ToolBar'
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
} from 'components/FlowChart/DialogContext'
import { FileSelectDialog } from 'components/common/FileSelectDialog'
import { FormHelperText, Popover } from '@mui/material'
import { NodeData } from 'store/slice/FlowElement/FlowElementType'

const initDialogFile = {
  filePath: '',
  open: false,
  fileTreeType: undefined,
  multiSelect: false,
  onSelectFile: () => null,
}

const ReactFlowProviderComponent = ReactFlowProvider as FC<{
  children: ReactNode
}>

export const ReactFlowComponent = React.memo<UseRunPipelineReturnType>(
  (props) => {
    const flowElements = useSelector(selectFlowElements)
    const nodes = flowElements.filter((item) =>
      isNode(item),
    ) as Node<NodeData>[]
    const egdes = flowElements.filter(
      (item) => !isNode(item),
    ) as Edge<NodeData>[]
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
        setFlowElements([
          ...nodes,
          ...addEdge(
            {
              ...params,
              animated: false,
              style: { width: 5 },
              type: 'buttonedge',
            },
            egdes,
          ),
        ]),
      )
    }

    const onNodesChange = (changes: NodeChange[]) => {
      dispatch(setFlowElements([...applyNodeChanges(changes, nodes), ...egdes]))
    }

    const onDragOver = (event: DragEvent) => {
      event.preventDefault()
      event.dataTransfer.dropEffect = 'move'
    }

    const onNodeDragStop = (event: ReactMouseEvent, node: Node) => {
      dispatch(
        editFlowElementPositionById({
          nodeId: node.id,
          coord: { x: node.position.x, y: node.position.y },
        }),
      )
    }

    const flowPosition = useSelector(selectFlowPosition)

    const onMoveEnd = (_: MouseEvent | TouchEvent, viewport: Viewport) => {
      dispatch(setFlowPosition(viewport))
    }

    const [reactFlowInstance, setReactFlowInstance] = React.useState<
      ReactFlowInstance<NodeData, NodeData> | undefined
    >()

    const onLoad = (
      reactFlowInstance: ReactFlowInstance<NodeData, NodeData>,
    ) => {
      setReactFlowInstance(reactFlowInstance)
    }

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
          <ReactFlowProviderComponent>
            <div className="reactflow-wrapper" ref={wrapparRef}>
              <ReactFlow
                ref={drop}
                nodes={nodes}
                edges={egdes}
                onNodesChange={onNodesChange}
                onConnect={onConnect}
                onInit={onLoad}
                onDragOver={onDragOver}
                onNodeDragStop={onNodeDragStop}
                nodeTypes={reactFlowNodeTypes as unknown as NodeTypes}
                edgeTypes={reactFlowEdgeTypes as unknown as EdgeTypes}
                defaultViewport={flowPosition}
                onMoveEnd={onMoveEnd}
              >
                <ToolBar {...props} />
                <Controls />
              </ReactFlow>
            </div>
          </ReactFlowProviderComponent>
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
