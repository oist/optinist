import { memo } from "react"
import { useSelector } from "react-redux"
import { useNavigate } from "react-router-dom"

import ArrowBackIosIcon from "@mui/icons-material/ArrowBackIos"
import { Button } from "@mui/material"
import Box from "@mui/material/Box"

import { CreateWorkflowButton } from "components/Workspace/FlowChart/Buttons/ClearWorkflow"
import { ImportWorkflowConfigButton } from "components/Workspace/FlowChart/Buttons/ImportWorkflowConfig"
import { NWBSettingButton } from "components/Workspace/FlowChart/Buttons/NWB"
import { RunButtons } from "components/Workspace/FlowChart/Buttons/RunButtons"
import { SnakemakeButton } from "components/Workspace/FlowChart/Buttons/Snakemake"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"
import { selectModeStandalone } from "store/slice/Standalone/StandaloneSeclector"

export const ToolBar = memo(function ToolBar(props: UseRunPipelineReturnType) {
  const navigate = useNavigate()
  const isStandalone = useSelector(selectModeStandalone)
  return (
    <Box
      style={{
        display: "flex",
        alignItems: "center",
        position: "absolute",
        float: "right",
        textAlign: "right",
        top: -7,
        right: 10,
        zIndex: 4,
        textTransform: "none",
        fontSize: "1rem",
      }}
    >
      {!isStandalone && (
        <Button onClick={() => navigate("/console/workspaces")}>
          <ArrowBackIosIcon />
          Workspaces
        </Button>
      )}
      <CreateWorkflowButton />
      <ImportWorkflowConfigButton />
      <SnakemakeButton />
      <NWBSettingButton />
      <RunButtons {...props} />
    </Box>
  )
})
