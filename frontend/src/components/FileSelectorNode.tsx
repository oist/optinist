import React, { CSSProperties } from 'react'
import { useDispatch } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { uploadImageFile } from 'redux/slice/ImageIndex/ImageIndexAction'

export const FileSelectorNode = React.memo<NodeProps>((element) => {
  const targetHandleStyle: CSSProperties = { background: '#555' }
  const sourceHandleStyle: CSSProperties = { ...targetHandleStyle }
  const dispatch = useDispatch()

  const onFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    event.preventDefault()
    if (event.target.files != null && event.target.files[0] != null) {
      const file = event.target.files[0]
      const formData = new FormData()
      formData.append('file', file)
      const elementId = element.id
      const fileName = file.name
      dispatch(uploadImageFile({ elementId, fileName, formData }))
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
})
