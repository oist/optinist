import React, { memo, FC, CSSProperties } from 'react'

import { Handle, Position, NodeProps } from 'react-flow-renderer'

const ColorSelectorNode: FC<NodeProps> = () => {
  const targetHandleStyle: CSSProperties = { background: '#555' }
  const sourceHandleStyle: CSSProperties = { ...targetHandleStyle }

  // const inputRef = useRef<HTMLInputElement>(null);

  const [selectedFileName, setSelectedFileName] = React.useState<string | null>(
    null,
  )
  const onFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    if (event.target.files != null && event.target.files[0] != null) {
      console.log(event.target.files[0].name)
      setSelectedFileName(event.target.files[0].name)
    }
  }

  return (
    <>
      <input type="file" onChange={onFileInputChange} />
      {!!selectedFileName && selectedFileName}
      <Handle
        type="source"
        position={Position.Bottom}
        id="a"
        style={sourceHandleStyle}
      />
    </>
  )
}

export default memo(ColorSelectorNode)
