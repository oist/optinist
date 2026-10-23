import { memo, useState } from "react"
import { useDispatch, useSelector } from "react-redux"
import { Handle, Position, NodeProps } from "reactflow"

import { Action, ThunkAction } from "@reduxjs/toolkit"

import { StructureItemSelectDialog } from "components/Workspace/FlowChart/Dialog/StructureItemSelectDialog"
import { FileSelect } from "components/Workspace/FlowChart/FlowChartNode/FileSelect"
import { toHandleId } from "components/Workspace/FlowChart/FlowChartNode/FlowChartUtils"
import { NodeContainer } from "components/Workspace/FlowChart/FlowChartNode/NodeContainer"
import { HANDLE_STYLE } from "const/flowchart"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import { selectInputNodeDefined } from "store/slice/InputNode/InputNodeSelectors"
import { RootState } from "store/store"

export type TreeNodeType = TreeDirType | TreeFileType

export interface TreeDirType {
  path: string
  name: string
  isDir: true
  nodes: TreeNodeType[]
  dataType?: string | null
}

export interface TreeFileType {
  path: string
  name: string
  isDir: false
  dataType?: string | null
  shape?: number[] | null
  nbytes?: string
}

export interface FileNodeConfig {
  fileType: string
  handleId: string
  handleType: string
  treeKeyPrefix: string
  selectFilePath: (
    nodeId: string,
  ) => (state: RootState) => string | string[] | undefined
  selectStructurePath: (
    nodeId: string,
  ) => (state: RootState) => string | undefined
  setStructurePath: (params: {
    nodeId: string
    path: string
  }) => Action<unknown>
  getTree: (params: {
    path: string
    workspaceId: number
  }) => ThunkAction<unknown, RootState, unknown, Action<unknown>>
  selectTree: () => (state: RootState) => TreeNodeType[] | undefined
  selectIsLoading: () => (state: RootState) => boolean
}

export function createStructuredFileNode(config: FileNodeConfig) {
  const FileNode = memo(function FileNode(element: NodeProps) {
    const defined = useSelector(selectInputNodeDefined(element.id))
    if (defined) {
      return <FileNodeImple {...element} config={config} />
    } else {
      return null
    }
  })
  FileNode.displayName = `${config.fileType}FileNode`
  return FileNode
}

const FileNodeImple = memo(function FileNodeImple({
  id: nodeId,
  selected,
  config,
}: NodeProps & { config: FileNodeConfig }) {
  const dispatch = useDispatch()
  const filePathRaw = useSelector(config.selectFilePath(nodeId))
  const filePath = Array.isArray(filePathRaw) ? filePathRaw[0] : filePathRaw

  const [open, setOpen] = useState(false)
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowNodeById(nodeId))
  }

  return (
    <NodeContainer nodeId={nodeId} selected={selected}>
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: "black", position: "absolute", top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        nodeId={nodeId}
        onChangeFilePath={(path) => {
          if (!Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        setOpen={setOpen}
        fileType={config.fileType}
        filePath={filePath ?? ""}
      />
      {filePath !== undefined && (
        <StructureItemSelectDialog
          open={open}
          setOpen={setOpen}
          nodeId={nodeId}
          config={config}
        />
      )}
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, config.handleId, config.handleType)}
        style={{ ...HANDLE_STYLE }}
      />
    </NodeContainer>
  )
})
