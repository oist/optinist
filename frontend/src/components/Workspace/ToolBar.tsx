import { memo } from "react"

import Box from "@mui/material/Box"

import { CreateWorkflowButton } from "components/Workspace/FlowChart/Buttons/CreateWorkflow"
import { ImportWorkflowConfigButton } from "components/Workspace/FlowChart/Buttons/ImportWorkflowConfig"
import { NWBSettingButton } from "components/Workspace/FlowChart/Buttons/NWB"
import { RunButtons } from "components/Workspace/FlowChart/Buttons/RunButtons"
import { SnakemakeButton } from "components/Workspace/FlowChart/Buttons/Snakemake"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"

export const ToolBar = memo(function ToolBar(props: UseRunPipelineReturnType) {
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
      <CreateWorkflowButton />
      <ImportWorkflowConfigButton />
      <SnakemakeButton />
      <NWBSettingButton />
      <RunButtons {...props} />
    </Box>
  )
})
