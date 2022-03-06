import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@mui/material/Button'
import Box from '@mui/material/Box'
import Alert from '@mui/material/Alert'
import AlertTitle from '@mui/material/AlertTitle'
import IconButton from '@mui/material/IconButton'
import Table from '@mui/material/Table'
import TableBody from '@mui/material/TableBody'
import TableCell, { tableCellClasses } from '@mui/material/TableCell'
import TableContainer from '@mui/material/TableContainer'
import TableHead from '@mui/material/TableHead'
import TableRow from '@mui/material/TableRow'
import Paper from '@mui/material/Paper'
import Dialog from '@mui/material/Dialog'
import DialogActions from '@mui/material/DialogActions'
import DialogTitle from '@mui/material/DialogTitle'
import KeyboardArrowDownIcon from '@mui/icons-material/KeyboardArrowDown'
import KeyboardArrowUpIcon from '@mui/icons-material/KeyboardArrowUp'
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline'
import GetAppIcon from '@mui/icons-material/GetApp'
import ReplayIcon from '@mui/icons-material/Replay'
import { useSnackbar } from 'notistack'

import { CollapsibleTable } from './CollapsibleTable'
import {
  selectExperimentsSatusIsUninitialized,
  selectExperimentsSatusIsFulfilled,
  selectExperimentUidList,
  selectExperimentTimeStamp,
  selectExperimentName,
  selectExperimentStatus,
  selectExperimentsSatusIsError,
  selectExperimentsErrorMessage,
} from 'store/slice/Experiments/ExperimentsSelectors'
import {
  deleteExperimentByUid,
  importExperimentByUid,
} from 'store/slice/Experiments/ExperimentsActions'
import { getExperiments } from 'store/slice/Experiments/ExperimentsActions'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { ExperimentStatusIcon } from './ExperimentStatusIcon'
import { AppDispatch } from 'store/store'

export const ExperimentTable: React.FC = () => {
  const isUninitialized = useSelector(selectExperimentsSatusIsUninitialized)
  const isFulfilled = useSelector(selectExperimentsSatusIsFulfilled)
  const isError = useSelector(selectExperimentsSatusIsError)
  const dispatch = useDispatch()
  React.useEffect(() => {
    if (isUninitialized) {
      dispatch(getExperiments())
    }
  }, [dispatch, isUninitialized])
  if (isFulfilled) {
    return <TableImple />
  } else if (isError) {
    return <ExperimentsErrorView />
  } else {
    return null
  }
}

const ExperimentsErrorView: React.FC = () => {
  const message = useSelector(selectExperimentsErrorMessage)
  return (
    <Alert severity="error">
      <AlertTitle>failed to get experiments...</AlertTitle>
      {message}
    </Alert>
  )
}

export const ExperimentUidContext = React.createContext<string>('')

const TableImple: React.FC = () => {
  const uidList = useSelector(selectExperimentUidList, arrayEqualityFn)
  const dispatch = useDispatch()
  const onClickReload = () => {
    dispatch(getExperiments())
  }
  return (
    <Box>
      <Box
        sx={{
          display: 'flex',
          justifyContent: 'flex-end',
        }}
      >
        <Button
          sx={{
            marginBottom: (theme) => theme.spacing(1),
          }}
          variant="outlined"
          endIcon={<ReplayIcon />}
          onClick={onClickReload}
        >
          Reload
        </Button>
      </Box>
      <TableContainer component={Paper} elevation={0} variant="outlined">
        <Table aria-label="collapsible table">
          <TableHead>
            <Head />
          </TableHead>
          <TableBody>
            {uidList.map((uid) => (
              <ExperimentUidContext.Provider value={uid} key={uid}>
                <Row key={uid} />
              </ExperimentUidContext.Provider>
            ))}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  )
}

const Head: React.FC = () => {
  return (
    <TableRow>
      <TableCell />
      <TableCell>Timestamp</TableCell>
      <TableCell>Name</TableCell>
      <TableCell>Success</TableCell>
      <TableCell>Import</TableCell>
      <TableCell>Delete</TableCell>
    </TableRow>
  )
}

const Row: React.FC = () => {
  const uid = React.useContext(ExperimentUidContext)
  const timestamp = useSelector(selectExperimentTimeStamp(uid))
  const status = useSelector(selectExperimentStatus(uid))
  const name = useSelector(selectExperimentName(uid))
  const [open, setOpen] = React.useState(false)
  return (
    <React.Fragment>
      <TableRow
        sx={{
          '& > *': {
            borderBottom: 'unset',
          },
          [`& .${tableCellClasses.root}`]: {
            borderBottomWidth: 0,
          },
        }}
      >
        <TableCell>
          <IconButton
            aria-label="expand row"
            size="small"
            onClick={() => setOpen((prevOpen) => !prevOpen)}
          >
            {open ? <KeyboardArrowUpIcon /> : <KeyboardArrowDownIcon />}
          </IconButton>
        </TableCell>
        <TableCell component="th" scope="row">
          {timestamp}
        </TableCell>
        <TableCell>{name}</TableCell>
        <TableCell>
          <ExperimentStatusIcon status={status} />
        </TableCell>
        <TableCell>
          <ImportExperimentButton />
        </TableCell>
        <TableCell>
          <DeleteExperimentButton />
        </TableCell>
      </TableRow>
      <CollapsibleTable open={open} />
    </React.Fragment>
  )
}

const ImportExperimentButton: React.FC = () => {
  const dispatch: AppDispatch = useDispatch()
  const uid = React.useContext(ExperimentUidContext)
  const { enqueueSnackbar } = useSnackbar()

  const onClick = () => {
    dispatch(importExperimentByUid(uid))
      .unwrap()
      .then(() =>
        enqueueSnackbar('Successfully imported.', { variant: 'success' }),
      )
  }
  return (
    <IconButton onClick={onClick}>
      <GetAppIcon color="primary" />
    </IconButton>
  )
}

const DeleteExperimentButton: React.FC = () => {
  const dispatch = useDispatch()
  const uid = React.useContext(ExperimentUidContext)

  const name = useSelector(selectExperimentName(uid))
  const [open, setOpen] = React.useState(false)

  const onClickOpen = () => {
    setOpen(true)
  }
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    setOpen(false)
    dispatch(deleteExperimentByUid(uid))
  }
  return (
    <>
      <IconButton onClick={onClickOpen}>
        <DeleteOutlineIcon color="error" />
      </IconButton>
      <Dialog open={open}>
        <DialogTitle>Are you sure you want to delete {name}?</DialogTitle>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined" color="inherit">
            Cancel
          </Button>
          <Button onClick={onClickOk} variant="outlined" autoFocus>
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
}
