import React, { ChangeEvent, useState } from 'react'
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
import TablePagination from '@mui/material/TablePagination'
import Paper from '@mui/material/Paper'
import KeyboardArrowDownIcon from '@mui/icons-material/KeyboardArrowDown'
import KeyboardArrowUpIcon from '@mui/icons-material/KeyboardArrowUp'
import ReplayIcon from '@mui/icons-material/Replay'
import DeleteIcon from '@mui/icons-material/Delete'
import Checkbox from '@mui/material/Checkbox'
import Dialog from '@mui/material/Dialog'
import DialogActions from '@mui/material/DialogActions'
import DialogTitle from '@mui/material/DialogTitle'
import TableSortLabel from '@mui/material/TableSortLabel'
import Typography from '@mui/material/Typography'

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
  selectExperimentHasNWB,
} from 'store/slice/Experiments/ExperimentsSelectors'
import {
  deleteExperimentByList,
  getExperiments,
} from 'store/slice/Experiments/ExperimentsActions'
import { ExperimentStatusIcon } from './ExperimentStatusIcon'
import { Experiment } from 'store/slice/Experiments/ExperimentsType'
import { DeleteButton } from './Button/DeleteButton'
import {
  NWBDownloadButton,
  ConfigDownloadButton,
} from './Button/DownloadButton'
import { ImportButton } from './Button/ImportButton'
import { useLocalStorage } from 'components/utils/LocalStorageUtil'
import { styled } from '@mui/material/styles'
import { renameExperiment } from 'api/experiments/Experiments'

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

const LOCAL_STORAGE_KEY_PER_PAGE = 'optinist_experiment_table_per_page'

const TableImple = React.memo(() => {
  const experimentList = useSelector(selectExperimentList)
  const experimentListValues = Object.values(experimentList)
  const experimentListKeys = Object.keys(experimentList)
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

  const onChangeAllCheck = (checked: boolean) => {
    if (checked) {
      setCheckedList(experimentListValues.map((experiment) => experiment.uid))
    } else {
      setCheckedList([])
    }
  }

  const recordsIsEmpty = experimentListKeys.length === 0

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

  const [page, setPage] = React.useState(0)

  const handleChangePage = (event: unknown, newPage: number) => {
    setPage(newPage)
  }

  const [rowsPerPage, setRowsPerPage] = useLocalStorage(
    LOCAL_STORAGE_KEY_PER_PAGE,
    10,
    (value) => {
      const valueNum = Number(value)
      return isNaN(valueNum) ? 10 : valueNum
    },
  )
  const handleChangeRowsPerPage = (
    event: React.ChangeEvent<HTMLInputElement>,
  ) => {
    const newValue = parseInt(event.target.value, 10)
    setRowsPerPage(newValue)
    setPage(0)
  }

  // Avoid a layout jump when reaching the last page with empty rows.
  const emptyRows =
    page > 0
      ? Math.max(0, (1 + page) * rowsPerPage - experimentListKeys.length)
      : 0

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column' }}>
      <Box
        sx={{
          display: 'flex',
          justifyContent: 'flex-end',
          alignItems: 'center',
        }}
      >
        {!recordsIsEmpty && (
          <Typography sx={{ flexGrow: 1, m: 1 }}>
            {checkedList.length} selected
          </Typography>
        )}
        <Button
          sx={{
            margin: (theme) => theme.spacing(0, 1, 1, 0),
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
      <Paper
        elevation={0}
        variant="outlined"
        sx={{
          flexGlow: 1,
          height: '100%',
        }}
      >
        <TableContainer component={Paper} elevation={0}>
          <Table aria-label="collapsible table">
            <HeadItem
              order={order}
              sortHandler={sortHandler}
              allCheckIndeterminate={
                checkedList.length !== 0 &&
                checkedList.length !== Object.keys(experimentList).length
              }
              allChecked={
                checkedList.length === Object.keys(experimentList).length
              }
              onChangeAllCheck={onChangeAllCheck}
              checkboxVisible={!recordsIsEmpty}
            />
            <TableBody>
              {experimentListValues
                .sort(getComparator(order, sortTarget))
                .slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage)
                .map((expData) => (
                  <ExperimentUidContext.Provider
                    value={expData.uid}
                    key={expData.uid}
                  >
                    <RowItem
                      onCheckBoxClick={onCheckBoxClick}
                      checked={checkedList.includes(expData.uid)}
                    />
                  </ExperimentUidContext.Provider>
                ))}
              {emptyRows > 0 && (
                <TableRow
                  style={{
                    height: 53 * emptyRows,
                  }}
                >
                  <TableCell colSpan={10} />
                </TableRow>
              )}
              {recordsIsEmpty && (
                <TableRow>
                  <TableCell colSpan={10}>
                    <Typography
                      sx={{
                        color: (theme) => theme.palette.text.secondary,
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        height: '300px',
                        textAlign: 'center',
                      }}
                      variant="h6"
                    >
                      No Rows...
                    </Typography>
                  </TableCell>
                </TableRow>
              )}
            </TableBody>
          </Table>
        </TableContainer>
        <TablePagination
          rowsPerPageOptions={[5, 10, 25]}
          component="div"
          count={experimentListKeys.length}
          rowsPerPage={rowsPerPage}
          page={page}
          onPageChange={handleChangePage}
          onRowsPerPageChange={handleChangeRowsPerPage}
        />
      </Paper>
    </Box>
  )
})

const HeadItem = React.memo<{
  order: Order
  sortHandler: any
  allChecked: boolean
  onChangeAllCheck: (checked: boolean) => void
  allCheckIndeterminate: boolean
  checkboxVisible: boolean
}>(
  ({
    order,
    sortHandler,
    allChecked,
    onChangeAllCheck,
    allCheckIndeterminate,
    checkboxVisible,
  }) => {
    return (
      <TableHead>
        <TableRow>
          <TableCell padding="checkbox">
            <Checkbox
              sx={{ visibility: checkboxVisible ? 'visible' : 'hidden' }}
              checked={allChecked}
              indeterminate={allCheckIndeterminate}
              onChange={(e) => onChangeAllCheck(e.target.checked)}
            />
          </TableCell>
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
            <TableSortLabel
              active
              direction={order}
              onClick={sortHandler('uid')}
            >
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
          <TableCell>SnakeFile</TableCell>
          <TableCell>NWB</TableCell>
          <TableCell>Delete</TableCell>
        </TableRow>
      </TableHead>
    )
  },
)

const RowItem = React.memo<{
  onCheckBoxClick: (uid: string) => void
  checked: boolean
}>(({ onCheckBoxClick, checked }) => {
  const uid = React.useContext(ExperimentUidContext)
  const timestamp = useSelector(selectExperimentTimeStamp(uid))
  const status = useSelector(selectExperimentStatus(uid))
  const name = useSelector(selectExperimentName(uid))
  const hasNWB = useSelector(selectExperimentHasNWB(uid))
  const [open, setOpen] = React.useState(false)
  const [isEdit, setEdit] = useState(false)
  const [error, setErrorEdit] = useState('')
  const [valueEdit, setValueEdit] = useState(name)
  const dispatch = useDispatch()

  const onBlur = (event: any) => {
    event.preventDefault()
    if (error) return
    setTimeout(() => {
      setEdit(false)
      onSaveNewName()
    }, 300)
  }

  const onEdit = (event: any) => {
    if (isEdit || error) return
    event.preventDefault()
    setEdit(true)
  }

  const onChangeName = (event: ChangeEvent<HTMLInputElement>) => {
    let errorEdit = ''
    if (!event.target.value.trim()) {
      errorEdit = 'Name is empty'
    }
    setErrorEdit(errorEdit)
    setValueEdit(event.target.value)
  }

  const onSaveNewName = async () => {
    if (valueEdit === name) return
    await renameExperiment(uid, valueEdit)
    dispatch(getExperiments())
  }

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
        <TableCell padding="checkbox">
          <Checkbox onChange={() => onCheckBoxClick(uid)} checked={checked} />
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
        <TableCell sx={{ width: 160, position: 'relative' }} onClick={onEdit}>
          {!isEdit ? (
            valueEdit
          ) : (
            <>
              <Input
                placeholder="Name"
                error={!!error}
                onChange={onChangeName}
                autoFocus
                onBlur={onBlur}
                value={valueEdit}
              />
              {error ? <TextError>{error}</TextError> : null}
            </>
          )}
        </TableCell>
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
          <NWBDownloadButton name={uid} hasNWB={hasNWB} />
        </TableCell>
        <TableCell>
          <DeleteButton />
        </TableCell>
      </TableRow>
      <CollapsibleTable open={open} />
    </React.Fragment>
  )
})

const Input = styled('input')<{ error: boolean }>(({ error }) => ({
  width: '100%',
  border: 'none',
  borderBottom: '1px solid',
  outline: 'none',
  color: error ? '#d32f2f' : '',
  borderColor: error ? '#d32f2f' : '',
}))

const TextError = styled(Typography)(() => ({
  color: '#d32f2f',
  fontSize: 12,
  height: 12,
  position: 'absolute',
  bottom: 12,
}))

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
