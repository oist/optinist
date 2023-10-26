import { memo } from "react"
import { useDispatch } from "react-redux"

import RouteIcon from "@mui/icons-material/Route"
import { IconButton, Tooltip } from "@mui/material"

import { toggleSnakemake } from "store/slice/RightDrawer/RightDrawerSlice"

export const SnakemakeButton = memo(function SnakemakeButton() {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleSnakemake())
  }
  return (
    <Tooltip title="Snakemake settings">
      <IconButton onClick={handleClick}>
        <RouteIcon color="primary" />
      </IconButton>
    </Tooltip>
  )
})
