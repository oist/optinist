import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { selectPipelineLatestUid } from 'store/slice/Pipeline/PipelineSelectors'
import { List, ListItem, ListSubheader, Divider } from '@mui/material'
import {
  selectExperimentName,
  selectExperimentsSatusIsFulfilled,
  selectExperimentsSatusIsUninitialized,
} from 'store/slice/Experiments/ExperimentsSelectors'
import { getExperiments } from 'store/slice/Experiments/ExperimentsActions'

export const CurrentPipelineInfo: React.FC = () => {
  const uid = useSelector(selectPipelineLatestUid)
  const isUninitialized = useSelector(selectExperimentsSatusIsUninitialized)
  const isFulFilled = useSelector(selectExperimentsSatusIsFulfilled)
  const dispatch = useDispatch()
  React.useEffect(() => {
    if (isUninitialized) {
      dispatch(getExperiments())
    }
  }, [dispatch, isUninitialized])

  return (
    <>
      {uid && (
        <>
          <List dense>
            <ExperimentUidInfo uid={uid} />
            {isFulFilled && <ExperimentNameInfo uid={uid} />}
          </List>
          <Divider />
        </>
      )}
    </>
  )
}

const ExperimentUidInfo = React.memo<{ uid: string }>(({ uid }) => {
  return (
    <>
      <ListSubheader>ID</ListSubheader>
      <ListItem>{uid}</ListItem>
    </>
  )
})

const ExperimentNameInfo = React.memo<{ uid: string }>(({ uid }) => {
  const name = useSelector(selectExperimentName(uid))
  return (
    <>
      <ListSubheader>NAME</ListSubheader>
      <ListItem>{name}</ListItem>
    </>
  )
})
