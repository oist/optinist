import { memo } from "react"

import Box from "@mui/material/Box"
import LinearProgress, {
  LinearProgressProps,
} from "@mui/material/LinearProgress"
import Typography from "@mui/material/Typography"

interface LinearProgressWithLabelProps extends LinearProgressProps {
  value: number
}

export const LinearProgressWithLabel = memo(function LinearProgressWithLabel(
  props: LinearProgressWithLabelProps,
) {
  return (
    <Box display="flex" alignItems="center">
      <Box width="100%" mr={1}>
        <LinearProgress variant="determinate" {...props} />
      </Box>
      <Box minWidth={35}>
        <Typography variant="body2" color="textSecondary">{`${Math.round(
          props.value,
        )}%`}</Typography>
      </Box>
    </Box>
  )
})
