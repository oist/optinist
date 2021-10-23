import React, { CSSProperties } from 'react'
import { useDispatch } from 'react-redux'
import { alpha, useTheme } from '@material-ui/core'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { uploadImageFile } from 'redux/slice/ImageIndex/ImageIndexAction'
import { FlexLayoutModelContext } from 'App'
import { useTabAction } from 'FlexLayoutHook'
import { OUTPUT_TABSET_ID } from 'const/flexlayout'

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
  const model = React.useContext(FlexLayoutModelContext)
  const actionForImageTab = useTabAction(element.id, 'image', OUTPUT_TABSET_ID)
  const onClick = () => {
    if (actionForImageTab != null) {
      model.doAction(actionForImageTab)
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
      onClick={onClick}
    >
      <div style={{ padding: 5 }}>
        <input type="file" onChange={onFileInputChange} />
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
