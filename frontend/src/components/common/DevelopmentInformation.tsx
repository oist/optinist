import { FC, memo, useState, useEffect } from "react"

import { Typography, Grid } from "@mui/material"

import { ParamSection } from "components/common/ParamSection"

interface Versions {
  version: string
  nodeVersion: string
  reactVersion: string
}
export const DevelopmentInformation: FC = () => {
  const [versions, setVersions] = useState<Versions | null>(null)

  useEffect(() => {
    const fetchData = async () => {
      try {
        // Note:
        // In consideration of the case where .versions.json does not exist,
        // dynamic import is performed using require.
        //
        // eslint-disable-next-line @typescript-eslint/no-var-requires
        const versionsJson = require(".versions.json") as Versions

        setVersions(versionsJson)
      } catch (err) {
        // set empty values
        setVersions({ version: "-", nodeVersion: "-", reactVersion: "-" })
      }
    }

    fetchData()
  }, [])

  return (
    <>
      <ParamSection title="Development Info &#x2699;">
        <Grid container mb={2}>
          <LabelValueGridRow
            label="Studio ver."
            value={versions?.version || ""}
          />
          <LabelValueGridRow
            label="Node ver."
            value={versions?.nodeVersion || ""}
          />
          <LabelValueGridRow
            label="React ver."
            value={versions?.reactVersion || ""}
          />
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
