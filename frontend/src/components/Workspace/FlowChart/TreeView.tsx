import { memo, useCallback, useEffect } from "react"
import { useDrag } from "react-dnd"
import { useSelector, useDispatch } from "react-redux"

import AddIcon from "@mui/icons-material/Add"
import ChevronRightIcon from "@mui/icons-material/ChevronRight"
import ExpandMoreIcon from "@mui/icons-material/ExpandMore"
import { styled, Typography } from "@mui/material"
import IconButton from "@mui/material/IconButton"
import { treeItemClasses } from "@mui/x-tree-view"
import { TreeItem } from "@mui/x-tree-view/TreeItem"
import { TreeView } from "@mui/x-tree-view/TreeView"

import {
  DND_ITEM_TYPE_SET,
  TreeItemCollectedProps,
  TreeItemDragObject,
  TreeItemDropResult,
} from "components/Workspace/FlowChart/DnDItemType"
import { REACT_FLOW_NODE_TYPE, REACT_FLOW_NODE_TYPE_KEY } from "const/flowchart"
import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"
import {
  selectAlgorithmListIsLatest,
  selectAlgorithmListTree,
} from "store/slice/AlgorithmList/AlgorithmListSelectors"
import {
  AlgorithmChild,
  AlgorithmNodeType,
} from "store/slice/AlgorithmList/AlgorithmListType"
import { isAlgoChild } from "store/slice/AlgorithmList/AlgorithmListUtils"
import {
  addAlgorithmNode,
  addInputNode,
} from "store/slice/FlowElement/FlowElementActions"
import {
  NODE_TYPE,
  NODE_TYPE_SET,
} from "store/slice/FlowElement/FlowElementType"
import { FILE_TYPE, FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { selectPipelineLatestUid } from "store/slice/Pipeline/PipelineSelectors"
import { AppDispatch } from "store/store"
import { getNanoId } from "utils/nanoid/NanoIdUtils"

export const AlgorithmTreeView = memo(function AlgorithmTreeView() {
  const dispatch = useDispatch<AppDispatch>()
  const algoList = useSelector(selectAlgorithmListTree)
  const isLatest = useSelector(selectAlgorithmListIsLatest)
  const workflowId = useSelector(selectPipelineLatestUid)
  const runAlready = typeof workflowId !== "undefined"

  useEffect(() => {
    if (!isLatest) {
      dispatch(getAlgoList())
    }
  }, [dispatch, isLatest])

  const onAddAlgoNode = useCallback(
    (
      nodeName: string,
      functionPath: string,
      position?: { x: number; y: number },
    ) => {
      const name = nodeName
      const newNode = {
        id: `${name}_${getNanoId()}`,
        type: REACT_FLOW_NODE_TYPE_KEY.AlgorithmNode,
        data: { label: name, type: NODE_TYPE_SET.ALGORITHM },
        position,
      }
      dispatch(
        addAlgorithmNode({
          node: newNode,
          name,
          functionPath,
          runAlready,
        }),
      )
    },
    [dispatch, runAlready],
  )

  return (
    <TreeView
      sx={{
        flexGrow: 1,
        height: "100%",
      }}
      defaultCollapseIcon={<ExpandMoreIcon />}
      defaultExpandIcon={<ChevronRightIcon />}
    >
      <TreeItem nodeId="Data" label="Data">
        <InputNodeComponent
          fileName={"image"}
          nodeName={"imageData"}
          fileType={FILE_TYPE_SET.IMAGE}
        />
        <InputNodeComponent
          fileName={"csv"}
          nodeName={"csvData"}
          fileType={FILE_TYPE_SET.CSV}
        />
        <InputNodeComponent
          fileName={"hdf5"}
          nodeName={"hdf5Data"}
          fileType={FILE_TYPE_SET.HDF5}
        />
        <InputNodeComponent
          fileName={"fluo"}
          nodeName={"fluoData"}
          fileType={FILE_TYPE_SET.FLUO}
        />
        <InputNodeComponent
          fileName={"behavior"}
          nodeName={"behaviorData"}
          fileType={FILE_TYPE_SET.BEHAVIOR}
        />
        <InputNodeComponent
          fileName={"matlab"}
          nodeName={"matlabData"}
          fileType={FILE_TYPE_SET.MATLAB}
        />
        <InputNodeComponent
          fileName={"microscope"}
          nodeName={"microscopeData"}
          fileType={FILE_TYPE_SET.MICROSCOPE}
        />
      </TreeItem>
      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.entries(algoList).map(([name, node], i) => (
          <AlgoNodeComponentRecursive
            name={name}
            node={node}
            onAddAlgoNode={onAddAlgoNode}
            key={i.toFixed()}
          />
        ))}
      </TreeItem>
    </TreeView>
  )
})

interface InputNodeComponentProps {
  fileName: string
  nodeName: string
  fileType: FILE_TYPE
}

const InputNodeComponent = memo(function InputNodeComponent({
  fileName,
  nodeName,
  fileType,
}: InputNodeComponentProps) {
  const dispatch = useDispatch()

  const onAddDataNode = useCallback(
    (
      nodeType: NODE_TYPE,
      nodeName: string,
      fileType: FILE_TYPE,
      position?: { x: number; y: number },
    ) => {
      let reactFlowNodeType: REACT_FLOW_NODE_TYPE | "" = ""
      switch (fileType) {
        case FILE_TYPE_SET.CSV:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.CsvFileNode
          break
        case FILE_TYPE_SET.IMAGE:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.ImageFileNode
          fileType = FILE_TYPE_SET.IMAGE
          break
        case FILE_TYPE_SET.HDF5:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.HDF5FileNode
          fileType = FILE_TYPE_SET.HDF5
          break
        case FILE_TYPE_SET.FLUO:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.FluoFileNode
          fileType = FILE_TYPE_SET.FLUO
          break
        case FILE_TYPE_SET.BEHAVIOR:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.BehaviorFileNode
          fileType = FILE_TYPE_SET.BEHAVIOR
          break
        case FILE_TYPE_SET.MATLAB:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.MatlabFileNode
          fileType = FILE_TYPE_SET.MATLAB
          break
        case FILE_TYPE_SET.MICROSCOPE:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.MicroscopeFileNode
          fileType = FILE_TYPE_SET.MICROSCOPE
          break
      }
      const newNode = {
        id: `input_${getNanoId()}`,
        type: reactFlowNodeType,
        data: { label: nodeName, type: nodeType },
        position,
      }
      dispatch(addInputNode({ node: newNode, fileType }))
    },
    [dispatch],
  )

  const { isDragging, dragRef } = useLeafItemDrag(
    useCallback(
      (position) => {
        onAddDataNode(NODE_TYPE_SET.INPUT, nodeName, fileType, position)
      },
      [onAddDataNode, nodeName, fileType],
    ),
  )

  return (
    <LeafItem
      ref={dragRef}
      style={{
        opacity: isDragging ? 0.6 : 1,
      }}
      onFocusCapture={(e) => e.stopPropagation()}
      nodeId={fileName}
      label={
        <AddButton
          name={fileName}
          onClick={() => onAddDataNode(NODE_TYPE_SET.INPUT, nodeName, fileType)}
        />
      }
    />
  )
})

interface AlgoNodeComponentBaseProps {
  name: string
  onAddAlgoNode: (
    nodeName: string,
    functionPath: string,
    position?: { x: number; y: number },
  ) => void
}

interface AlgoNodeComponentRecursiveProps extends AlgoNodeComponentBaseProps {
  node: AlgorithmNodeType
}

const AlgoNodeComponentRecursive = memo(function AlgoNodeComponentRecursive({
  name,
  node,
  onAddAlgoNode,
}: AlgoNodeComponentRecursiveProps) {
  if (isAlgoChild(node)) {
    return (
      <AlgoNodeComponent
        name={name}
        node={node}
        onAddAlgoNode={onAddAlgoNode}
      />
    )
  } else {
    return (
      <TreeItem nodeId={name} label={name}>
        {Object.entries(node.children).map(([name, node], i) => (
          <AlgoNodeComponentRecursive
            name={name}
            node={node}
            key={i.toFixed()}
            onAddAlgoNode={onAddAlgoNode}
          />
        ))}
      </TreeItem>
    )
  }
})

interface AlgoNodeComponentProps extends AlgoNodeComponentBaseProps {
  node: AlgorithmChild
}

const AlgoNodeComponent = memo(function AlgoNodeComponent({
  name,
  node,
  onAddAlgoNode,
}: AlgoNodeComponentProps) {
  const { isDragging, dragRef } = useLeafItemDrag(
    useCallback(
      (position) => {
        onAddAlgoNode(name, node.functionPath, position)
      },
      [onAddAlgoNode, name, node],
    ),
  )
  return (
    <LeafItem
      ref={dragRef}
      style={{
        opacity: isDragging ? 0.6 : 1,
      }}
      onFocusCapture={(e) => e.stopPropagation()}
      nodeId={name}
      label={
        <AddButton
          name={name}
          onClick={() => onAddAlgoNode(name, node.functionPath)}
        />
      }
    />
  )
})

interface AddButtonProps {
  name: string
  onClick: () => void
}

const AddButton = memo(function AddButton({ name, onClick }: AddButtonProps) {
  return (
    <>
      <IconButton
        aria-label="add"
        style={{ padding: 2 }}
        size="large"
        onClick={onClick}
      >
        <AddIcon />
      </IconButton>
      <Typography
        variant="inherit"
        style={{
          textOverflow: "ellipsis",
          overflow: "visible",
          width: "8rem",
          display: "inline-block",
        }}
      >
        {name}
      </Typography>
    </>
  )
})

// 未使用icon分の幅を消す
const LeafItem = styled(TreeItem)({
  // background: 'red',
  [`& .${treeItemClasses.iconContainer}`]: {
    margin: 0,
    width: 0,
  },
})

function useLeafItemDrag(
  onDragEnd: (position: { x: number; y: number }) => void,
) {
  const [{ isDragging }, dragRef] = useDrag<
    TreeItemDragObject,
    TreeItemDropResult,
    TreeItemCollectedProps
  >(
    () => ({
      type: DND_ITEM_TYPE_SET.TREE_ITEM,
      end: (_, monitor) => {
        const position = monitor.getDropResult()?.position
        if (monitor.didDrop() && position != null) {
          onDragEnd(position)
        }
      },
      collect: (monitor) => ({
        isDragging: monitor.isDragging(),
      }),
    }),
    [onDragEnd],
  )
  return { isDragging, dragRef }
}
