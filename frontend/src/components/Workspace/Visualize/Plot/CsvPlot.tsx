import { memo, useContext, useEffect, useMemo } from "react"
import { useSelector, useDispatch } from "react-redux"

import { LinearProgress, Typography } from "@mui/material"
import { DataGrid, GridColDef } from "@mui/x-data-grid"

import { CsvData } from "api/outputs/Outputs"
import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getCsvData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectCsvData,
  selectCsvDataError,
  selectCsvDataIsFulfilled,
  selectCsvDataIsInitialized,
  selectCsvDataIsPending,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectCsvItemSetHeader,
  selectCsvItemSetIndex,
  selectCsvItemTranspose,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"
import { twoDimarrayEqualityFn } from "utils/EqualityUtils"

export const CsvPlot = memo(function CsvPlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const isInitialized = useSelector(selectCsvDataIsInitialized(path))
  const isPending = useSelector(selectCsvDataIsPending(path))
  const isFulfilled = useSelector(selectCsvDataIsFulfilled(path))
  const error = useSelector(selectCsvDataError(path))
  const dispatch = useDispatch<AppDispatch>()
  useEffect(() => {
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

const CsvPlotImple = memo(function CsvPlotImple() {
  const { itemId, filePath: path } = useContext(DisplayDataContext)
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
interface PresentationalCsvPlotProps {
  path: string
  transpose: boolean
  setHeader: number | null
  setIndex: boolean
}

export const PresentationalCsvPlot = memo(function PresentationalCsvPlot({
  path,
  transpose,
  setIndex,
  setHeader,
}: PresentationalCsvPlotProps) {
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

  const csvData = useMemo(() => {
    if (transpose) {
      return data[0].map((col, i) => data.map((row) => row[i]))
    }
    return data
  }, [data, transpose])

  const header: GridColDef[] = useMemo(() => {
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
        { field: "col0", headerName: "index", width: 150 },
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
      const rowObj = Object.fromEntries(
        [row_id, ...row].map((value, index) => {
          return [`col${index}`, value]
        }),
      )
      rowObj["id"] = row_id
      return rowObj
    })
    .filter(
      (value, idx) =>
        setHeader === null || (setHeader !== null && idx > setHeader),
    )

  return (
    <div style={{ height: 300, width: "100%" }}>
      <DataGrid rows={rows} columns={header} />
    </div>
  )
})
