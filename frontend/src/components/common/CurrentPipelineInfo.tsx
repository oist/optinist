import { FC, memo } from "react"
import { useSelector } from "react-redux"

import { Divider, Typography, Grid } from "@mui/material"

import {
  selectExperimentName,
  selectExperimentsStatusIsFulfilled,
} from "store/slice/Experiments/ExperimentsSelectors"
import { selectPipelineLatestUid } from "store/slice/Pipeline/PipelineSelectors"

export const CurrentPipelineInfo: FC = () => {
  const uid = useSelector(selectPipelineLatestUid)
  const isFulFilled = useSelector(selectExperimentsStatusIsFulfilled)

  return (
    <>
      {uid && (
        <>
          <Grid container paddingX={2} paddingBottom={1}>
            <ExperimentUidInfo uid={uid} />
            {isFulFilled && <ExperimentNameInfo uid={uid} />}
          </Grid>
          <Divider />
        </>
      )}
    </>
  )
}

interface UidProps {
  uid: string
}

const ExperimentUidInfo = memo(function ExperimentUidInfo({ uid }: UidProps) {
  return <LabelValueGridRow label="ID" value={uid} />
})

const ExperimentNameInfo = memo(function ExperimentNameInfo({ uid }: UidProps) {
  const name = useSelector(selectExperimentName(uid))
  return <LabelValueGridRow label="NAME" value={name} />
})

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
      <Grid item xs={4}>
        <Typography variant="body2" color="textSecondary">
          {label}:
        </Typography>
      </Grid>
      <Grid item xs={8}>
        <Typography variant="body2" color="textSecondary">
          {value}
        </Typography>
      </Grid>
    </>
  )
})
