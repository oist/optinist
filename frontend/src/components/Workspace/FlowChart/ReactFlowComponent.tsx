import {
  memo,
  DragEvent,
  FC,
  MouseEvent,
  ReactNode,
  useState,
  useRef,
} from "react"
import { useDrop } from "react-dnd"
import { useSelector, useDispatch } from "react-redux"
import {
  ReactFlow,
  isNode,
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
} from "reactflow"
import "reactflow/dist/style.css"

import "style/flow.css"
import { FormHelperText, Popover } from "@mui/material"

import { FileSelectDialog } from "components/common/FileSelectDialog"
import {
  DialogContext,
  ErrorDialogValue,
  OpenDialogValue,
} from "components/Workspace/FlowChart/DialogContext"
import {
  DND_ITEM_TYPE_SET,
  TreeItemCollectedProps,
  TreeItemDragObject,
  TreeItemDropResult,
} from "components/Workspace/FlowChart/DnDItemType"
import { AlgorithmOutputDialog } from "components/Workspace/FlowChart/FlowChartNode/AlgorithmOutputDialog"
import {
  reactFlowEdgeTypes,
  reactFlowNodeTypes,
} from "components/Workspace/FlowChart/FlowChartNode/ReactFlowNodeTypesConst"
import { ToolBar } from "components/Workspace/ToolBar"
import {
  selectFlowEdges,
  selectFlowNodes,
  selectFlowPosition,
} from "store/slice/FlowElement/FlowElementSelectors"
import {
  editFlowNodePositionById,
  setEdgesChange,
  setFlowEdges,
  setFlowPosition,
  setNodesChange,
} from "store/slice/FlowElement/FlowElementSlice"
import { NodeData } from "store/slice/FlowElement/FlowElementType"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"

const initDialogFile = {
  filePath: "",
  open: false,
  fileTreeType: undefined,
  multiSelect: false,
  onSelectFile: () => null,
}

const ReactFlowProviderComponent = ReactFlowProvider as FC<{
  children: ReactNode
}>

export const ReactFlowComponent = memo(function ReactFlowComponent(
  props: UseRunPipelineReturnType,
) {
  const flowNodes = useSelector(selectFlowNodes)
  const flowEdges = useSelector(selectFlowEdges)
  const egdes = flowEdges.filter((item) => !isNode(item)) as Edge<NodeData>[]
  const dispatch = useDispatch()
  const [dialogNodeId, setDialogNodeId] = useState("")
  const [dialogFile, setDialogFile] = useState<OpenDialogValue>(initDialogFile)
  const [messageError, setMessageError] = useState<ErrorDialogValue>({
    anchorElRef: { current: null },
    message: "",
  })

  const onConnect = (params: Connection | Edge) => {
    dispatch(
      setFlowEdges([
        ...addEdge(
          {
            ...params,
            animated: false,
            style: { width: 5 },
            type: "buttonedge",
          },
          egdes,
        ),
      ]),
    )
  }

  const onNodesChange: OnNodesChange = (changes) =>
    dispatch(setNodesChange(changes))
  const onEdgesChange: OnEdgesChange = (changes) =>
    dispatch(setEdgesChange(changes))

  const onDragOver = (event: DragEvent) => {
    event.preventDefault()
    event.dataTransfer.dropEffect = "move"
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
    useState<ReactFlowInstance>()

  const onInit = (reactFlowInstance: ReactFlowInstance) =>
    setReactFlowInstance(reactFlowInstance)

  const wrapparRef = useRef<HTMLDivElement>(null)
  const [, drop] = useDrop<
    TreeItemDragObject,
    TreeItemDropResult,
    TreeItemCollectedProps
  >(
    () => ({
      accept: DND_ITEM_TYPE_SET.TREE_ITEM,
      drop: (_, monitor) => {
        let position: TreeItemDropResult["position"] = undefined
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
              defaultViewport={{
                x: flowPosition[0],
                y: flowPosition[1],
                zoom: flowPosition[2],
              }}
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
            onClose={() => setDialogNodeId("")}
          />
        )}
        {dialogFile.open && (
          <FileSelectDialog
            multiSelect={dialogFile.multiSelect}
            initialFilePath={dialogFile.filePath}
            open={dialogFile.open}
            onClickOk={(path: string | string[]) => {
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
                message: "",
              })
            }
            anchorOrigin={{
              vertical: "top",
              horizontal: "right",
            }}
            transformOrigin={{
              vertical: "bottom",
              horizontal: "left",
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
})
