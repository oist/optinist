import React from "react"
import { useNavigate } from "react-router-dom"

import ArrowBackIosIcon from "@mui/icons-material/ArrowBackIos"
import { Button } from "@mui/material"
import Box from "@mui/material/Box"

import { ClearWorkflowButton } from "components/Workspace/FlowChart/Buttons/ClearWorkflow"
import { ImportWorkflowConfigButton } from "components/Workspace/FlowChart/Buttons/ImportWorkflowConfig"
import { NWBSettingButton } from "components/Workspace/FlowChart/Buttons/NWB"
import { RunButtons } from "components/Workspace/FlowChart/Buttons/RunButtons"
import { SnakemakeButton } from "components/Workspace/FlowChart/Buttons/Snakemake"
import { IS_STANDALONE } from "const/Mode"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"




export const ToolBar = React.memo<UseRunPipelineReturnType>((props) => {
  const navigate = useNavigate()
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
      {!IS_STANDALONE && (
        <Button onClick={() => navigate("/console/workspaces")}>
          <ArrowBackIosIcon />
          Workspaces
        </Button>
      )}
      <ImportWorkflowConfigButton />
      <ClearWorkflowButton />
      <SnakemakeButton />
      <NWBSettingButton />
      <RunButtons {...props} />
    </Box>
  )
})
