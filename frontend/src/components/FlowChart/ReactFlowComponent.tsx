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
  FlowTransform,
} from 'react-flow-renderer'

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
import { ImageFileNode } from './FlowChartNode/ImageFileNode'
import { AlgorithmNode } from './FlowChartNode/AlgorithmNode'
import { CsvFileNode } from './FlowChartNode/CsvFileNode'
import { HDF5FileNode } from './FlowChartNode/HDF5FileNode'
import { CustomEdge } from './CustomEdge'

const componentTypes = {
  ImageFileNode,
  CsvFileNode,
  HDF5FileNode,
  AlgorithmNode,
} as const

const edgeTypes = {
  buttonedge: CustomEdge,
} as const

export const ReactFlowComponent = React.memo(() => {
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

  return (
    <div className="flow">
      <ReactFlowProvider>
        <div className="reactflow-wrapper">
          <ReactFlow
            elements={flowElements}
            onElementsRemove={onElementsRemove}
            onConnect={onConnect}
            onDragOver={onDragOver}
            onNodeDragStop={onNodeDragStop}
            nodeTypes={componentTypes}
            edgeTypes={edgeTypes}
            defaultPosition={[flowPosition.x, flowPosition.y]}
            defaultZoom={flowPosition.zoom}
            onMoveEnd={onMoveEnd}
          >
            <Controls />
          </ReactFlow>
        </div>
      </ReactFlowProvider>
    </div>
  )
})
