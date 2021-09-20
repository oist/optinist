import { memo, FC, CSSProperties } from 'react'
import { useDispatch } from 'react-redux'
import { setElementPath } from 'redux/slice/Element/Element'
import { Handle, Position, NodeProps } from 'react-flow-renderer'

const FileSelectorNode: FC<NodeProps> = (element) => {
  const targetHandleStyle: CSSProperties = { background: '#555' }
  const sourceHandleStyle: CSSProperties = { ...targetHandleStyle }
  const dispatch = useDispatch()

  // const dispatch = useDispatch()

  const onFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    if (event.target.files != null && event.target.files[0] != null) {
      // console.log(event.target.files[0].name)
      // dispatch(setSelectedFileName(event.target.files[0].name));
      dispatch(
        setElementPath({ id: element.id, path: event.target.files[0].name }),
      )
    }
  }

  return (
    <>
      <input type="file" onChange={onFileInputChange} />
      <Handle
        type="source"
        position={Position.Bottom}
        id="a"
        style={sourceHandleStyle}
      />
    </>
  )
}

export default memo(FileSelectorNode)
