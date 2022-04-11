import React, { useState } from 'react'
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
import KeyboardArrowDownIcon from '@mui/icons-material/KeyboardArrowDown'
import KeyboardArrowUpIcon from '@mui/icons-material/KeyboardArrowUp'
import ReplayIcon from '@mui/icons-material/Replay'
import DeleteIcon from '@mui/icons-material/Delete'

import { CollapsibleTable } from './CollapsibleTable'
import {
  selectExperimentsSatusIsUninitialized,
  selectExperimentsSatusIsFulfilled,
  selectExperimentTimeStamp,
  selectExperimentName,
  selectExperimentStatus,
  selectExperimentsSatusIsError,
  selectExperimentsErrorMessage,
  selectExperimentList,
} from 'store/slice/Experiments/ExperimentsSelectors'
import {
  deleteExperimentByList,
  getExperiments,
} from 'store/slice/Experiments/ExperimentsActions'
import { ExperimentStatusIcon } from './ExperimentStatusIcon'
import { Experiment } from 'store/slice/Experiments/ExperimentsType'
import {
  Checkbox,
  Dialog,
  DialogActions,
  DialogTitle,
  TableSortLabel,
} from '@mui/material'
import { DeleteButton } from './Button/DeleteButton'
import {
  NWBDownloadButton,
  ConfigDownloadButton,
} from './Button/DownloadButton'
import { ImportButton } from './Button/ImportButton'

export const ExperimentUidContext = React.createContext<string>('')

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

const TableImple = React.memo(() => {
  const experimentList = useSelector(selectExperimentList)
  const dispatch = useDispatch()
  const onClickReload = () => {
    dispatch(getExperiments())
  }
  const [order, setOrder] = React.useState<Order>('asc')
  const [sortTarget, setSortTarget] =
    React.useState<keyof Experiment>('timestamp')
  const sortHandler =
    (property: keyof Experiment) => (event: React.MouseEvent<unknown>) => {
      const isAsc = sortTarget === property && order === 'asc'
      setOrder(isAsc ? 'desc' : 'asc')
      setSortTarget(property)
    }

  const [checkedList, setCheckedList] = useState<string[]>([])
  const [open, setOpen] = React.useState(false)

  const onCheckBoxClick = (uid: string) => {
    if (checkedList.includes(uid)) {
      setCheckedList(checkedList.filter((v) => v !== uid))
    } else {
      setCheckedList([...checkedList, uid])
    }
  }
  const onClickDelete = () => {
    setOpen(true)
  }
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    dispatch(deleteExperimentByList(checkedList))
    setCheckedList([])
    setOpen(false)
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
        <Button
          sx={{
            marginBottom: (theme) => theme.spacing(1),
          }}
          variant="outlined"
          color="error"
          endIcon={<DeleteIcon />}
          onClick={onClickDelete}
          disabled={checkedList.length === 0}
        >
          Delete
        </Button>
      </Box>
      <Dialog open={open}>
        <DialogTitle>Are you sure you want to delete?</DialogTitle>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined" color="inherit">
            Cancel
          </Button>
          <Button onClick={onClickOk} variant="outlined" autoFocus>
            OK
          </Button>
        </DialogActions>
      </Dialog>
      <TableContainer component={Paper} elevation={0} variant="outlined">
        <Table aria-label="collapsible table">
          <HeadItem order={order} sortHandler={sortHandler} />
          <TableBody>
            {Object.values(experimentList)
              .sort(getComparator(order, sortTarget))
              .map((expData) => (
                <ExperimentUidContext.Provider
                  value={expData.uid}
                  key={expData.uid}
                >
                  <RowItem onCheckBoxClick={onCheckBoxClick} />
                </ExperimentUidContext.Provider>
              ))}
          </TableBody>
        </Table>
      </TableContainer>
    </Box>
  )
})

const HeadItem = React.memo<{
  order: Order
  sortHandler: any
}>(({ order, sortHandler }) => {
  return (
    <TableHead>
      <TableRow>
        <TableCell />
        <TableCell />
        <TableCell>
          <TableSortLabel
            active
            direction={order}
            onClick={sortHandler('timestamp')}
          >
            Timestamp
          </TableSortLabel>
        </TableCell>
        <TableCell>
          <TableSortLabel active direction={order} onClick={sortHandler('uid')}>
            ID
          </TableSortLabel>
        </TableCell>
        <TableCell>
          <TableSortLabel
            active
            direction={order}
            onClick={sortHandler('name')}
          >
            Name
          </TableSortLabel>
        </TableCell>
        <TableCell>Success</TableCell>
        <TableCell>Reproduce</TableCell>
        <TableCell>workflow</TableCell>
        <TableCell>NWB</TableCell>
        <TableCell>Delete</TableCell>
      </TableRow>
    </TableHead>
  )
})

const RowItem = React.memo<{
  onCheckBoxClick: (uid: string) => void
}>(({ onCheckBoxClick }) => {
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
          <Checkbox onChange={() => onCheckBoxClick(uid)} />
        </TableCell>
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
        <TableCell>{uid}</TableCell>
        <TableCell>{name}</TableCell>
        <TableCell>
          <ExperimentStatusIcon status={status} />
        </TableCell>
        <TableCell>
          <ImportButton />
        </TableCell>
        <TableCell>
          <ConfigDownloadButton />
        </TableCell>
        <TableCell>
          <NWBDownloadButton name={uid} />
        </TableCell>
        <TableCell>
          <DeleteButton />
        </TableCell>
      </TableRow>
      <CollapsibleTable open={open} />
    </React.Fragment>
  )
})

type Order = 'asc' | 'desc'

function getComparator<Key extends keyof any>(
  order: Order,
  orderBy: Key,
): (
  a: { [key in Key]: number | string },
  b: { [key in Key]: number | string },
) => number {
  return order === 'desc'
    ? (a, b) => descendingComparator(a, b, orderBy)
    : (a, b) => -descendingComparator(a, b, orderBy)
}

function descendingComparator<T>(a: T, b: T, orderBy: keyof T) {
  if (b[orderBy] < a[orderBy]) {
    return -1
  }
  if (b[orderBy] > a[orderBy]) {
    return 1
  }
  return 0
}
