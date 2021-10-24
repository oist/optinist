import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { uploadImageFile } from 'redux/slice/ImageIndex/ImageIndexAction'
import { alpha, useTheme } from '@material-ui/core'
import { runStatusSelector } from 'redux/slice/Element/ElementSelector'

export const FileSelectorNode = React.memo<NodeProps>((element) => {
  const targetHandleStyle: CSSProperties = {
    width: 8,
    height: '100%',
    border: '1px solid',
    borderColor: '#555',
    borderRadius: 0,
  }
  const sourceHandleStyle: CSSProperties = { ...targetHandleStyle }
  const dispatch = useDispatch()

  const runStatus = useSelector(runStatusSelector)
  const inputRef = React.useRef<HTMLInputElement>(null)
  const [filePathError, setFilePathError] = React.useState(false)

  React.useEffect(() => {
    if (inputRef.current != null) {
      if (runStatus === 'failed' && inputRef.current.files?.length === 0) {
        setFilePathError(true)
      }
    }
  }, [runStatus, inputRef.current])

  const onFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    event.preventDefault()
    if (event.target.files != null && event.target.files[0] != null) {
      const file = event.target.files[0]
      const formData = new FormData()
      formData.append('file', file)
      const elementId = element.id
      const fileName = file.name
      dispatch(uploadImageFile({ elementId, fileName, formData }))
      setFilePathError(false)
    }
  }
  const theme = useTheme()
  return (
    <div
      style={{
        height: '100%',
        background: element.selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
    >
      <div
        style={{
          padding: 5,
          color: filePathError ? theme.palette.error.main : undefined,
        }}
      >
        <input ref={inputRef} type="file" onChange={onFileInputChange} />
      </div>
      <Handle
        type="source"
        position={Position.Right}
        id="a"
        style={sourceHandleStyle}
      />
    </div>
  )
})
