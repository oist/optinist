import React from 'react'
import { useSelector } from 'react-redux'
import { selectPipelineLatestUid } from 'store/slice/Pipeline/PipelineSelectors'
import { Divider, Typography, Grid } from '@mui/material'
import {
  selectExperimentName,
  selectExperimentsSatusIsFulfilled,
} from 'store/slice/Experiments/ExperimentsSelectors'

export const CurrentPipelineInfo: React.FC = () => {
  const uid = useSelector(selectPipelineLatestUid)
  const isFulFilled = useSelector(selectExperimentsSatusIsFulfilled)

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

const ExperimentUidInfo = React.memo<{ uid: string }>(({ uid }) => {
  return <LabelValueGridRow label="ID" value={uid} />
})

const ExperimentNameInfo = React.memo<{ uid: string }>(({ uid }) => {
  const name = useSelector(selectExperimentName(uid))
  return <LabelValueGridRow label="NAME" value={name} />
})

const LabelValueGridRow = React.memo<{ label: string; value: string }>(
  ({ label, value }) => {
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
  },
)
