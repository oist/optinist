import { useDispatch } from 'react-redux'
import { EdgeProps, getBezierPath } from 'reactflow'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import 'style/flowbutton.css'

const foreignObjectSize = 40

export const CustomEdge: React.FC<EdgeProps> = ({
  id,
  sourceX,
  sourceY,
  targetX,
  targetY,
  sourcePosition,
  targetPosition,
  style = {},
  markerEnd,
}) => {
  const [path, edgeCenterX, edgeCenterY] = getBezierPath({
    sourceX,
    sourceY,
    sourcePosition,
    targetX,
    targetY,
    targetPosition,
  })

  const dispatch = useDispatch()

  const onEdgeClick = () => {
    dispatch(deleteFlowElementsById(id))
  }

  return (
    <>
      <path
        id={id}
        style={style}
        className="react-flow__edge-path"
        d={path}
        markerEnd={markerEnd}
      />
      <foreignObject
        width={foreignObjectSize}
        height={foreignObjectSize}
        x={edgeCenterX - foreignObjectSize / 2}
        y={edgeCenterY - foreignObjectSize / 2}
        className="edgebutton-foreignobject"
      >
        <button className="flowbutton" onClick={onEdgeClick}>
          ×
        </button>
      </foreignObject>
    </>
  )
}
