import { FC, memo } from "react"

import { Typography, Grid } from "@mui/material"

import versions from ".versions.json"
import { ParamSection } from "components/common/ParamSection"

export const DevelopmentInformation: FC = () => {
  return (
    <>
      <ParamSection title="Development Info &#x2699;">
        <Grid container mb={2}>
          <LabelValueGridRow label="Studio ver." value={versions.version} />
          <LabelValueGridRow label="Node ver." value={versions.nodeVersion} />
          <LabelValueGridRow label="React ver." value={versions.reactVersion} />
        </Grid>
      </ParamSection>
    </>
  )
}

interface LabelValueGridRowProps {
  label: string
  value: string
}

const LabelValueGridRow = memo(function LabelValueGridRow({
  label,
  value,
}: LabelValueGridRowProps) {
  return (
    <>
      <Grid item xs={6}>
        <Typography fontWeight="bold" fontSize="0.95rem" color="text.secondary">
          {label}
        </Typography>
      </Grid>
      <Grid item xs={6}>
        <div>
          <Typography
            aria-haspopup="true"
            noWrap
            variant="body2"
            overflow="hidden"
            textOverflow="ellipsis"
          >
            {value}
          </Typography>
        </div>
      </Grid>
    </>
  )
})
