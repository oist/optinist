import { FC, memo, MouseEvent, useEffect, useRef, useState } from "react"
import { useSelector } from "react-redux"

import { Typography, Grid, Popover } from "@mui/material"

import { ParamSection } from "components/common/ParamSection"
import {
  selectCurrentPipelineName,
  selectPipelineLatestUid,
} from "store/slice/Pipeline/PipelineSelectors"
import { selectModeStandalone } from "store/slice/Standalone/StandaloneSeclector"
import {
  selectCurrentWorkspaceId,
  selectCurrentWorkspaceName,
} from "store/slice/Workspace/WorkspaceSelector"

export const CurrentPipelineInfo: FC = () => {
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const workspaceName = useSelector(selectCurrentWorkspaceName)
  const workflowId = useSelector(selectPipelineLatestUid)
  const workflowName = useSelector(selectCurrentPipelineName)
  const isStandalone = useSelector(selectModeStandalone)

  return (
    <>
      {!isStandalone && (
        <ParamSection title="Workspace">
          <Grid container mb={1}>
            {workspaceId && (
              <LabelValueGridRow
                label="ID"
                value={workspaceId.toFixed()}
                section="workspace"
              />
            )}
            {workspaceName && (
              <LabelValueGridRow
                label="NAME"
                value={workspaceName}
                section="workspace"
              />
            )}
          </Grid>
        </ParamSection>
      )}
      {workflowId && (
        <ParamSection title="Workflow">
          <Grid container mb={1}>
            <LabelValueGridRow
              label="ID"
              value={workflowId}
              section="workflow"
            />
            {workflowName && (
              <LabelValueGridRow
                label="NAME"
                value={workflowName}
                section="workflow"
              />
            )}
          </Grid>
        </ParamSection>
      )}
    </>
  )
}

interface LabelValueGridRowProps {
  label: string
  value: string
  section: string
}

const LabelValueGridRow = memo(function LabelValueGridRow({
  label,
  value,
  section,
}: LabelValueGridRowProps) {
  const [isOverflow, setIsOverflow] = useState(false)
  const ref = useRef<HTMLElement>(null)

  useEffect(() => {
    const element = ref.current
    if (element) {
      setIsOverflow(element.scrollWidth > element.offsetWidth)
    }
  }, [ref])

  const [anchorEl, setAnchorEl] = useState<HTMLButtonElement | null>(null)
  const handlePopoverOpen = (event: MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget)
  }
  const handlePopoverClose = () => {
    setAnchorEl(null)
  }
  const open = Boolean(anchorEl)
  const id = open ? `pipeline-info-popover-${section}-${label}` : undefined

  return (
    <>
      <Grid item xs={3}>
        <Typography fontWeight="bold" fontSize="0.95rem" color="text.secondary">
          {label}
        </Typography>
      </Grid>
      <Grid item xs={9}>
        <div>
          <Typography
            ref={ref}
            aria-owns={id}
            aria-haspopup="true"
            onMouseEnter={handlePopoverOpen}
            onMouseLeave={handlePopoverClose}
            noWrap
            variant="body2"
            overflow="hidden"
            textOverflow="ellipsis"
          >
            {value}
          </Typography>
          {isOverflow && (
            <Popover
              id={id}
              open={open}
              anchorEl={anchorEl}
              anchorOrigin={{
                vertical: "bottom",
                horizontal: "left",
              }}
              onClose={handlePopoverClose}
              sx={{
                pointerEvents: "none",
              }}
              disableRestoreFocus
            >
              <Typography variant="body2" padding={2}>
                {value}
              </Typography>
            </Popover>
          )}
        </div>
      </Grid>
    </>
  )
})
