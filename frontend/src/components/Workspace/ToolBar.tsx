import React from "react"
import Box from "@mui/material/Box"

import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"
import { NWBSettingButton } from "./FlowChart/Buttons/NWB"
import { SnakemakeButton } from "./FlowChart/Buttons/Snakemake"
import { RunButtons } from "./FlowChart/Buttons/RunButtons"
import { Button } from "@mui/material"
import ArrowBackIosIcon from "@mui/icons-material/ArrowBackIos"
import { useNavigate } from "react-router-dom"
import { IS_STANDALONE } from "const/Mode"
import { ImportWorkflowConfigButton } from "./FlowChart/Buttons/ImportWorkflowConfig"
import { ClearWorkflowButton } from "./FlowChart/Buttons/ClearWorkflow"

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
