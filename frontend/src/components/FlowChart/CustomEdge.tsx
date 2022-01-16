import {
  EdgeProps,
  getBezierPath,
  getMarkerEnd,
  getEdgeCenter,
} from 'react-flow-renderer'
import CloseOutlinedIcon from '@material-ui/icons/CloseOutlined'
import IconButton from '@mui/material/IconButton'
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
  arrowHeadType,
  markerEndId,
}) => {
  const edgePath = getBezierPath({
    sourceX,
    sourceY,
    sourcePosition,
    targetX,
    targetY,
    targetPosition,
  })
  const markerEnd = getMarkerEnd(arrowHeadType, markerEndId)
  const [edgeCenterX, edgeCenterY] = getEdgeCenter({
    sourceX,
    sourceY,
    targetX,
    targetY,
  })

  const onEdgeClick = (event: any, id: string) => {
    event.stopPropagation()
    alert(`remove ${id}`)
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
        // requiredExtensions="http://www.w3.org/1999/xhtml"
      >
        <body>
          <IconButton aria-label="delete">
            <CloseOutlinedIcon onClick={(event) => onEdgeClick(event, id)} />
          </IconButton>
        </body>
      </foreignObject>
    </>
  )
}
