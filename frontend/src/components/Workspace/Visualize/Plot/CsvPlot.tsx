import React from 'react'
import { LinearProgress, Typography } from '@mui/material'
import { useSelector, useDispatch } from 'react-redux'

import { DisplayDataContext } from '../DataContext'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectCsvData,
  selectCsvDataError,
  selectCsvDataIsFulfilled,
  selectCsvDataIsInitialized,
  selectCsvDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getCsvData } from 'store/slice/DisplayData/DisplayDataActions'
import { CsvData } from 'api/outputs/Outputs'
import { DataGrid, GridColDef } from '@mui/x-data-grid'
import {
  selectCsvItemSetHeader,
  selectCsvItemSetIndex,
  selectCsvItemTranspose,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { selectCurrentWorkspaceId } from 'store/slice/Workspace/WorkspaceSelector'

export const CsvPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const isInitialized = useSelector(selectCsvDataIsInitialized(path))
  const isPending = useSelector(selectCsvDataIsPending(path))
  const isFulfilled = useSelector(selectCsvDataIsFulfilled(path))
  const error = useSelector(selectCsvDataError(path))
  const dispatch = useDispatch()
  React.useEffect(() => {
    if (workspaceId && !isInitialized) {
      dispatch(getCsvData({ path, workspaceId }))
    }
  }, [dispatch, isInitialized, path, workspaceId])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <CsvPlotImple />
  } else {
    return null
  }
})

const CsvPlotImple = React.memo(() => {
  const { itemId, filePath: path } = React.useContext(DisplayDataContext)
  const transpose = useSelector(selectCsvItemTranspose(itemId))
  const setHeader = useSelector(selectCsvItemSetHeader(itemId))
  const setIndex = useSelector(selectCsvItemSetIndex(itemId))
  return (
    <PresentationalCsvPlot
      path={path}
      transpose={transpose}
      setHeader={setHeader}
      setIndex={setIndex}
    />
  )
})

/**
 * CsvFileNodeの設定時に表示するプレビューでも使用することを想定
 *
 * DisplayDataContextに依存しない
 */
export const PresentationalCsvPlot = React.memo<{
  path: string
  transpose: boolean
  setHeader: number | null
  setIndex: boolean
}>(({ path, transpose, setIndex, setHeader }) => {
  const data = useSelector(
    selectCsvData(path),
    (a: CsvData | undefined, b: CsvData | undefined) => {
      if (a != null && b != null) {
        return twoDimarrayEqualityFn(a, b)
      } else {
        return a === undefined && b === undefined
      }
    },
  )

  const csvData = React.useMemo(() => {
    if (transpose) {
      return data[0].map((col, i) => data.map((row) => row[i]))
    }
    return data
  }, [data, transpose])

  const header: GridColDef[] = React.useMemo(() => {
    const headerData = () => {
      if (setHeader === null) {
        return csvData[0]
      } else {
        if (setHeader >= csvData.length) {
          return csvData[csvData.length - 1]
        } else {
          return csvData[setHeader]
        }
      }
    }

    if (setIndex) {
      return [
        { field: 'col0', headerName: 'index', width: 150 },
        ...headerData().map((value, idx) => {
          return {
            field: `col${idx + 1}`,
            headerName: `${setHeader === null ? idx : value}`,
            width: 150,
          }
        }),
      ]
    } else {
      return headerData().map((value, idx) => {
        return {
          field: `col${idx + 1}`,
          headerName: `${setHeader === null ? idx : value}`,
          width: 150,
        }
      })
    }
  }, [csvData, setHeader, setIndex])
  const rows = csvData
    .map((row, row_id) => {
      let rowObj = Object.fromEntries(
        [row_id, ...row].map((value, index) => {
          return [`col${index}`, value]
        }),
      )
      rowObj['id'] = row_id
      return rowObj
    })
    .filter(
      (value, idx) =>
        setHeader === null || (setHeader !== null && idx > setHeader),
    )

  return (
    <div style={{ height: 300, width: '100%' }}>
      <DataGrid rows={rows} columns={header} />
    </div>
  )
})
