import { useDispatch } from 'react-redux'
import { EdgeProps, getBezierPath, getEdgeCenter } from 'react-flow-renderer'
import { deleteFlowEdgeById } from 'store/slice/FlowElement/FlowElementSlice'
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
  data,
  markerEnd,
}) => {
  const edgePath = getBezierPath({
    sourceX,
    sourceY,
    sourcePosition,
    targetX,
    targetY,
    targetPosition,
  })

  const [edgeCenterX, edgeCenterY] = getEdgeCenter({
    sourceX,
    sourceY,
    targetX,
    targetY,
  })

  const dispatch = useDispatch()

  const onEdgeClick = () => {
    dispatch(deleteFlowEdgeById(id))
  }

  return (
    <>
      <path
        id={id}
        style={style}
        className="react-flow__edge-path"
        d={edgePath}
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
