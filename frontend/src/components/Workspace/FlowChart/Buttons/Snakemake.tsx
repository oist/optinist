import { memo } from "react"
import { useDispatch, useSelector } from "react-redux"

import RouteIcon from "@mui/icons-material/Route"
import { IconButton, Tooltip } from "@mui/material"

import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"
import { toggleSnakemake } from "store/slice/RightDrawer/RightDrawerSlice"

export const SnakemakeButton = memo(function SnakemakeButton() {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleSnakemake())
  }
  const isPending = useSelector(selectPipelineIsStartedSuccess)
  return (
    <Tooltip title="Snakemake settings">
      <IconButton
        onClick={handleClick}
        color={"primary"}
        disabled={!!isPending}
      >
        <RouteIcon />
      </IconButton>
    </Tooltip>
  )
})
