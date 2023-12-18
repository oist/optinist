import { memo, useContext, useEffect, useMemo } from "react"
import { useSelector, useDispatch } from "react-redux"

import { LinearProgress, Typography } from "@mui/material"
import { DataGrid, GridColDef } from "@mui/x-data-grid"

import { MatlabData } from "api/outputs/Outputs"
import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { getMatlabData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectMatlabData,
  selectMatlabDataError,
  selectMatlabDataIsFulfilled,
  selectMatlabDataIsInitialized,
  selectMatlabDataIsPending,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectMatlabItemSetHeader,
  selectMatlabItemSetIndex,
  selectMatlabItemTranspose,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"
import { twoDimarrayEqualityFn } from "utils/EqualityUtils"

export const MatlabPlot = memo(function MatlabPlot() {
  const { filePath: path } = useContext(DisplayDataContext)
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const isInitialized = useSelector(selectMatlabDataIsInitialized(path))
  const isPending = useSelector(selectMatlabDataIsPending(path))
  const isFulfilled = useSelector(selectMatlabDataIsFulfilled(path))
  const error = useSelector(selectMatlabDataError(path))
  const dispatch = useDispatch<AppDispatch>()
  useEffect(() => {
    if (workspaceId && !isInitialized) {
      dispatch(getMatlabData({ path, workspaceId }))
    }
  }, [dispatch, isInitialized, path, workspaceId])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <MatlabPlotImple />
  } else {
    return null
  }
})

const MatlabPlotImple = memo(function MatlabPlotImple() {
  const { itemId, filePath: path } = useContext(DisplayDataContext)
  const transpose = useSelector(selectMatlabItemTranspose(itemId))
  const setHeader = useSelector(selectMatlabItemSetHeader(itemId))
  const setIndex = useSelector(selectMatlabItemSetIndex(itemId))
  return (
    <PresentationalMatlabPlot
      path={path}
      transpose={transpose}
      setHeader={setHeader}
      setIndex={setIndex}
    />
  )
})

/**
 * MatlabFileNodeの設定時に表示するプレビューでも使用することを想定
 *
 * DisplayDataContextに依存しない
 */
interface PresentationalMatlabPlotProps {
  path: string
  transpose: boolean
  setHeader: number | null
  setIndex: boolean
}

export const PresentationalMatlabPlot = memo(function PresentationalMatlabPlot({
  path,
  transpose,
  setIndex,
  setHeader,
}: PresentationalMatlabPlotProps) {
  const data = useSelector(
    selectMatlabData(path),
    (a: MatlabData | undefined, b: MatlabData | undefined) => {
      if (a != null && b != null) {
        return twoDimarrayEqualityFn(a, b)
      } else {
        return a === undefined && b === undefined
      }
    },
  )

  const matlabData = useMemo(() => {
    if (transpose) {
      return data[0].map((col, i) => data.map((row) => row[i]))
    }
    return data
  }, [data, transpose])

  const header: GridColDef[] = useMemo(() => {
    const headerData = () => {
      if (setHeader === null) {
        return matlabData[0]
      } else {
        if (setHeader >= matlabData.length) {
          return matlabData[matlabData.length - 1]
        } else {
          return matlabData[setHeader]
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
  }, [matlabData, setHeader, setIndex])
  const rows = matlabData
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
