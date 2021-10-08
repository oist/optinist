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
      const uploadFolderName = `${file.name}(${element.id})`
      fetch(`http://localhost:8000/upload/${uploadFolderName}`, {
        method: 'POST',
        mode: 'cors',
        credentials: 'include',
        body: formData,
      })
        .then((response) => response.json())
        .then(
          (result) => {
            const elementId = element.id
            const fileName = file.name
            dispatch(
              uploadImageFile({
                elementId,
                fileName,
                folder: result.folderName,
                maxIndex: result.maxIndex,
              }),
            )
          },
          (error) => {
            console.log(error)
          },
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
})
