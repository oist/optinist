import React, { DragEvent, MouseEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import ReactFlow, {
  ReactFlowProvider,
  addEdge,
  Controls,
  Elements,
  Connection,
  Edge,
  Node,
  OnLoadParams,
  FlowTransform,
} from 'react-flow-renderer'
import { useDrop } from 'react-dnd'

import 'style/flow.css'
import {
  deleteFlowElements,
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

export const ReactFlowComponent = React.memo<UseRunPipelineReturnType>(
  (props) => {
    const flowElements = useSelector(selectFlowElements)
    const dispatch = useDispatch()

    const onConnect = (params: Connection | Edge) => {
      dispatch(
        setFlowElements(
          addEdge(
            {
              ...params,
              animated: false,
              style: { width: 5 },
              type: 'buttonedge',
            },
            flowElements,
          ),
        ),
      )
    }

    const onElementsRemove = (elementsToRemove: Elements) => {
      dispatch(deleteFlowElements(elementsToRemove))
    }

    const onDragOver = (event: DragEvent) => {
      event.preventDefault()
      event.dataTransfer.dropEffect = 'move'
    }

    const onNodeDragStop = (event: MouseEvent, node: Node) => {
      dispatch(
        editFlowElementPositionById({
          nodeId: node.id,
          coord: { x: node.position.x, y: node.position.y },
        }),
      )
    }

    const flowPosition = useSelector(selectFlowPosition)

    const onMoveEnd = (event: FlowTransform | undefined) => {
      if (event !== undefined) {
        dispatch(setFlowPosition(event))
      }
    }

    const [reactFlowInstance, setReactFlowInstance] =
      React.useState<OnLoadParams>()

    const onLoad = (reactFlowInstance: OnLoadParams) =>
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
        <ReactFlowProvider>
          <div className="reactflow-wrapper" ref={wrapparRef}>
            <ReactFlow
              ref={drop}
              elements={flowElements}
              onElementsRemove={onElementsRemove}
              onConnect={onConnect}
              onLoad={onLoad}
              onDragOver={onDragOver}
              onNodeDragStop={onNodeDragStop}
              nodeTypes={reactFlowNodeTypes}
              edgeTypes={reactFlowEdgeTypes}
              defaultPosition={[flowPosition.x, flowPosition.y]}
              defaultZoom={flowPosition.zoom}
              onMoveEnd={onMoveEnd}
            >
              <ToolBar {...props} />
              <Controls />
            </ReactFlow>
          </div>
        </ReactFlowProvider>
      </div>
    )
  },
)
