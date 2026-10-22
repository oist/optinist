import {
  memo,
  useContext,
  useCallback,
  useEffect,
  useState,
  MouseEvent,
  useRef,
  useMemo,
  ChangeEvent,
} from "react"
import PlotlyChart from "react-plotlyjs-ts"
import { useSelector, useDispatch } from "react-redux"

import createColormap from "colormap"
import { useSnackbar } from "notistack"
import {
  Datum,
  LayoutAxis,
  PlotData,
  PlotMouseEvent,
  PlotSelectionEvent,
  Annotations,
} from "plotly.js"

import { Button, LinearProgress, TextField, Typography } from "@mui/material"
import Box from "@mui/material/Box"
import FormControlLabel from "@mui/material/FormControlLabel"
import Slider from "@mui/material/Slider"
import { styled } from "@mui/material/styles"
import Switch from "@mui/material/Switch"

import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { useVisualize } from "components/Workspace/Visualize/VisualizeContext"
import {
  addRoi,
  cancelRoi,
  commitRoi,
  deleteRoi,
  mergeRoi,
  promoteRoi,
  getImageData,
  getRoiData,
  getStatus,
  getTimeSeriesInitData,
} from "store/slice/DisplayData/DisplayDataActions"
import {
  selectImageDataError,
  selectImageDataIsInitialized,
  selectImageDataIsPending,
  selectImageDataIsFulfilled,
  selectActiveImageData,
  selectRoiData,
  selectImageDataMaxSize,
  selectImageMeta,
  selectStatusRoi,
} from "store/slice/DisplayData/DisplayDataSelectors"
import {
  selectingImageArea,
  setImageItemClickedDataId,
} from "store/slice/VisualizeItem/VisualizeItemActions"
import {
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowLine,
  selectImageItemShowGrid,
  selectImageItemShowScale,
  selectImageItemColors,
  selectImageItemActiveIndex,
  selectImageItemStartIndex,
  selectImageItemEndIndex,
  selectRoiItemFilePath,
  selectRoiItemIndex,
  selectImageItemDuration,
  selectVisualizeItemWidth,
  selectVisualizeItemHeight,
  selectVisualizeSaveFilename,
  selectVisualizeSaveFormat,
  selectImageItemAlpha,
  selectVisualizeItems,
  selectImageItemShowRoiLabels,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  incrementImageActiveIndex,
  resetAllOrderList,
  setImageActiveIndex,
  setImageItemDuration,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import {
  isImageItem,
  isTimeSeriesItem,
} from "store/slice/VisualizeItem/VisualizeItemUtils"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { setDataCancel } from "store/slice/Workspace/WorkspaceSlice"
import { AppDispatch, RootState } from "store/store"
import { twoDimarrayEqualityFn } from "utils/EqualityUtils"

export type StatusROI = {
  temp_add_roi: number[]
  temp_delete_roi: number[]
  temp_merge_roi: number[]
  temp_promote_roi: number[]
}

const ADD_ROI = "Add ROI"
const DELETE_ROI = "Delete ROI"
const MERGE_ROI = "Merge ROI"
const PROMOTE_ROI = "Set as Cell ROI"
const CELL_ROI = "/cell_roi.json"
// all_roi is not offered: no commit path regenerates it
const PROMOTE_ROI_FILES = ["/noncell_roi.json", "/non_cell_roi.json"]

// the node types EditROI has commit wrappers for
const EDITABLE_ROI_NODES = ["suite2p", "lccd", "vacant_roi", "caiman"]

const isEditableNode = (path?: string | null) => {
  const nodeId = path?.split("/").slice(-2)[0] ?? ""
  return EDITABLE_ROI_NODES.some((name) => nodeId.includes(name))
}

const isEditableRoiPath = (path?: string | null): path is string =>
  isEditableNode(path) &&
  [CELL_ROI, ...PROMOTE_ROI_FILES].some((file) => path?.includes(file))
const WIDTH_CHARTJS = 321
const INIT_WIDTH_ROI = 30

const initSizeDrag = {
  width: INIT_WIDTH_ROI,
  height: INIT_WIDTH_ROI,
  left: Math.floor((WIDTH_CHARTJS - INIT_WIDTH_ROI) / 2),
  top: Math.floor((WIDTH_CHARTJS - INIT_WIDTH_ROI) / 2),
}

enum PositionDrag {
  "LEFT" = "LEFT",
  "RIGHT" = "RIGHT",
  "BOTTOM" = "BOTTOM",
  "TOP" = "TOP",
}

const sChart = 320

export const ImagePlot = memo(function ImagePlot() {
  const { filePath: path, itemId } = useContext(DisplayDataContext)

  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const startIndex = useSelector(selectImageItemStartIndex(itemId))
  const endIndex = useSelector(selectImageItemEndIndex(itemId))
  const isPending = useSelector(selectImageDataIsPending(path))
  const isInitialized = useSelector(selectImageDataIsInitialized(path))
  const isFulfilled = useSelector(selectImageDataIsFulfilled(path))
  const error = useSelector(selectImageDataError(path))

  const roiFilePath = useSelector(selectRoiItemFilePath(itemId))

  const dispatch = useDispatch<AppDispatch>()
  useEffect(() => {
    if (workspaceId) {
      if (!isInitialized) {
        dispatch(
          getImageData({
            path,
            workspaceId,
            startIndex: startIndex ?? 1,
            endIndex: endIndex ?? 10,
          }),
        )
      }
      if (roiFilePath != null) {
        dispatch(getRoiData({ path: roiFilePath, workspaceId }))
      }
    }
  }, [
    dispatch,
    isInitialized,
    path,
    workspaceId,
    startIndex,
    endIndex,
    roiFilePath,
  ])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <ImagePlotImple />
  } else {
    return null
  }
})

const ImagePlotImple = memo(function ImagePlotImple() {
  const { itemId } = useContext(DisplayDataContext)
  const activeIndex = useSelector(selectImageItemActiveIndex(itemId))
  return <ImagePlotChart activeIndex={activeIndex} />
})

interface ActiveIndexProps {
  activeIndex: number
}

const ImagePlotChart = memo(function ImagePlotChart({
  activeIndex,
}: ActiveIndexProps) {
  const { enqueueSnackbar } = useSnackbar()
  const dispatch = useDispatch<AppDispatch>()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const { filePath: path, itemId } = useContext(DisplayDataContext)
  const imageData = useSelector(
    selectActiveImageData(path, activeIndex),
    imageDataEqualtyFn,
  )
  const meta = useSelector(selectImageMeta(path))
  const roiFilePath = useSelector(selectRoiItemFilePath(itemId))

  const refRoiFilePath = useRef(roiFilePath)

  const roiData = useSelector(
    (state: RootState) =>
      roiFilePath != null ? selectRoiData(roiFilePath)(state) : [],
    imageDataEqualtyFn,
  )

  const [roiDataState, setRoiDataState] = useState(roiData)
  const itemsVisual = useSelector(selectVisualizeItems)
  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))
  const alpha = useSelector(selectImageItemAlpha(itemId))
  const timeDataMaxIndex = useSelector(selectRoiItemIndex(itemId, roiFilePath))
  const width = useSelector(selectVisualizeItemWidth(itemId))
  const height = useSelector(selectVisualizeItemHeight(itemId))
  const [sizeDrag, setSizeDrag] = useState(initSizeDrag)
  const [startDragAddRoi, setStartDragAddRoi] = useState(false)
  const [action, setAction] = useState("")
  const [positionDrag, setChangeSize] = useState<PositionDrag | undefined>()
  const showRoiLabels = useSelector(selectImageItemShowRoiLabels(itemId))

  const allowEditRoi = useMemo(
    () => isEditableNode(roiFilePath) && roiFilePath?.includes(CELL_ROI),
    [roiFilePath],
  )
  const allowPromoteRoi = useMemo(
    () =>
      isEditableNode(roiFilePath) &&
      PROMOTE_ROI_FILES.some((file) => roiFilePath?.includes(file)),
    [roiFilePath],
  )
  const { setRoisClick, roisClick, resetRoisClick } = useVisualize()
  const roiClicked = useMemo(() => roisClick[itemId] || [], [itemId, roisClick])

  const selectedStatus = useSelector(selectStatusRoi)

  const statusRoi = useMemo(() => {
    return (
      selectedStatus || {
        temp_add_roi: [],
        temp_delete_roi: [],
        temp_merge_roi: [],
        temp_promote_roi: [],
      }
    )
  }, [selectedStatus])

  const refPageXSize = useRef(0)
  const refPageYSize = useRef(0)

  function getRoiColor(roiIndex: number): number[] {
    const colors = createColormap({
      colormap: "jet",
      nshades: 200,
      format: "rgba",
    })
    return colors[(Math.abs(roiIndex) * 9) % 200]
  }

  useEffect(() => {
    setRoiDataState(roiData)
  }, [roiData])

  useEffect(() => {
    if (!roiFilePath || !workspaceId) return
    dispatch(getStatus({ path: roiFilePath, workspaceId }))
    //eslint-disable-next-line react-hooks/exhaustive-deps
  }, [roiFilePath, workspaceId])

  const refPromotePending = useRef(false)
  refPromotePending.current = statusRoi.temp_promote_roi.length > 0

  useEffect(() => {
    refRoiFilePath.current = roiFilePath
    return () => {
      const path = refRoiFilePath.current
      if (!workspaceId || !isEditableRoiPath(path)) return
      // a non-cell view cancels only when it has promotions of its own pending
      if (!path.includes(CELL_ROI) && !refPromotePending.current) return
      dispatch(cancelRoi({ path, workspaceId }))
    }
    //eslint-disable-next-line react-hooks/exhaustive-deps
  }, [roiFilePath])

  const data = useMemo(
    () => [
      {
        z: imageData,
        type: "heatmap",
        name: "images",
        colorscale: colorscale.map((value) => {
          let offset: number = parseFloat(value.offset)
          const offsets: number[] = colorscale.map((v) => {
            return parseFloat(v.offset)
          })
          // plotlyは端[0.0, 1.0]がないとダメなので、その設定
          if (offset === Math.max(...offsets)) {
            offset = 1.0
          }
          if (offset === Math.min(...offsets)) {
            offset = 0.0
          }
          const rgb = value.rgb
            .replace(/[^0-9,]/g, "")
            .split(",")
            .map((x) => Number(x))
          const hex = rgba2hex(rgb, alpha)
          return [offset, hex]
        }),
        // hoverinfo: isAddRoi || pointClick.length ? "none" : undefined,
        hoverongaps: false,
        showscale: showscale,
        zsmooth: zsmooth, // ["best", "fast", false]
      },
      {
        z: roiDataState,
        type: "heatmap",
        name: "roi",
        hovertemplate: action === ADD_ROI ? "none" : "ROI: %{z}",
        // hoverinfo: isAddRoi || pointClick.length ? "none" : undefined,
        colorscale: [...Array(timeDataMaxIndex + 1)].map((_, i) => {
          const offset: number = i / timeDataMaxIndex
          const rgba = getRoiColor(i)
          const hex = rgba2hex(rgba, 1)

          const isClickPoint = roiClicked.some((point) => point === i)
          const isDeleted = statusRoi?.temp_delete_roi?.includes(i) || false
          const isMerged = statusRoi?.temp_merge_roi?.includes(i) || false
          const isAdded = statusRoi?.temp_add_roi?.includes(i) || false
          const isPromoted = statusRoi?.temp_promote_roi?.includes(i) || false
          if ((allowEditRoi || allowPromoteRoi) && isClickPoint && action) {
            switch (action) {
              case DELETE_ROI:
                return [offset, "#FFA500"] // orange
              case MERGE_ROI:
                return [offset, "#e134eb"] // purple
              case ADD_ROI:
                return [offset, "3483eb"] // red
              case PROMOTE_ROI:
                return [offset, "#32A919"] // green
            }
          }
          if (isDeleted || isMerged || isAdded || isPromoted) {
            if (isDeleted) return [offset, "#FFA500"]
            if (isMerged) return [offset, "#e134eb"]
            if (isAdded) return [offset, "3483eb"]
            if (isPromoted) return [offset, "#32A919"]
          }
          if (isClickPoint || !roiClicked.length) return [offset, hex]
          return [offset, rgba2hex(rgba, 0.3)]
        }),
        zmin: 0,
        zmax: timeDataMaxIndex,
        hoverongaps: false,
        zsmooth: false,
        showscale: false,
      },
    ],
    [
      imageData,
      colorscale,
      showscale,
      zsmooth,
      roiDataState,
      timeDataMaxIndex,
      alpha,
      action,
      roiClicked,
      statusRoi,
      allowEditRoi,
      allowPromoteRoi,
    ],
  )

  const [selectMode, setSelectMode] = useState(false)
  const [cancelFirst, setCancelFirst] = useState(true)

  const checkStatus = useCallback(() => {
    if (!statusRoi) return
    return Object.keys(statusRoi).every(
      (key) => statusRoi[key as keyof StatusROI].length === 0,
    )
    //eslint-disable-next-line
  }, [JSON.stringify(statusRoi)])

  const [edit, setEdit] = useState<boolean>(!checkStatus())

  useEffect(() => {
    setEdit(!checkStatus())
  }, [checkStatus])

  const handleChange = (event: ChangeEvent<HTMLInputElement>) => {
    setSelectMode(event.target.checked)
  }

  const setRoiClickWhenSelect = (rois: string[]) => {
    resetRoisClick(itemId)
    rois.forEach((roi) => {
      setRoisClick(itemId, Number(roi))
      dispatch(
        setImageItemClickedDataId({
          itemId,
          clickedDataId: roi,
        }),
      )
    })
  }

  // debounceでイベントを間引きする。onSelectedはそれっぽい名前だが動かなかった。
  const onSelecting = debounce((event: PlotSelectionEvent) => {
    if (event.range != null) {
      dispatch(
        selectingImageArea({
          itemId,
          range: event.range,
          callback: setRoiClickWhenSelect,
        }),
      )
    }
  })

  const roiLabels = useMemo(() => {
    if (!showRoiLabels) return []
    const labels: Annotations[] = []
    const seen = new Set()

    roiDataState.forEach((row, _y) => {
      row.forEach((value, _x) => {
        if (value !== null && !seen.has(value)) {
          seen.add(value)
          const points = {
            x: [] as number[],
            y: [] as number[],
            count: 0,
          }
          roiDataState.forEach((r, yi) => {
            r.forEach((v, xi) => {
              if (v === value) {
                points.x.push(xi)
                points.y.push(yi)
                points.count++
              }
            })
          })

          const centerX = points.x.reduce((a, b) => a + b, 0) / points.count
          const centerY = points.y.reduce((a, b) => a + b, 0) / points.count

          const annotation: Partial<Annotations> = {
            x: centerX,
            y: centerY,
            text: `${value}`,
            xref: "x",
            yref: "y",
            showarrow: false,
            font: {
              color: "black", // Black text
              size: 10,
              weight: 700, // Bold text
            },
            bgcolor: "rgba(255, 255, 255, 0.6)", // Soft white semi-transparent background (60% opacity)
            borderpad: 0.5, // Padding around the text
          }

          labels.push(annotation as Annotations)
        }
      })
    })
    return labels
  }, [roiDataState, showRoiLabels])

  const layout = useMemo(
    () => ({
      title: {
        text: meta?.title,
        x: 0.1,
      },
      width: width,
      height: height - 130,
      margin: {
        t: 30, // top
        l: 100, // left
        b: 20, // bottom
      },
      dragmode: selectMode ? "select" : "pan",
      xaxis: {
        title: meta?.xlabel,
        autorange: true,
        showgrid: showgrid,
        showline: showline,
        zeroline: false,
        autotick: true,
        ticks: "",
        showticklabels: showticklabels,
        annotations: roiLabels,
      },
      yaxis: {
        title: meta?.ylabel,
        automargin: true,
        autorange: "reversed",
        showgrid: showgrid,
        showline: showline,
        zeroline: false,
        autotick: true,
        ticks: "",
        showticklabels: showticklabels,
      },
      annotations: roiLabels,
    }),
    //eslint-disable-next-line react-hooks/exhaustive-deps
    [
      meta,
      showgrid,
      showline,
      showticklabels,
      width,
      height,
      selectMode,
      action,
      roiLabels,
    ],
  )

  const saveFileName = useSelector(selectVisualizeSaveFilename(itemId))
  const saveFormat = useSelector(selectVisualizeSaveFormat(itemId))

  const config = {
    displayModeBar: true,
    responsive: true,
    toImageButtonOptions: {
      format: saveFormat,
      filename: saveFileName,
      // scale: number;
      // format: "png" | "svg" | "jpeg" | "webp";
      // height: number;
      // width: number;
    },
  }

  const onChartClick = (event: PlotMouseEvent) => {
    // use as unknown because original PlotDatum does not have z property
    const point: PlotDatum = event.points[0] as unknown as PlotDatum
    if (point.curveNumber >= 1 && point.z >= 0) {
      dispatch(
        setImageItemClickedDataId({ itemId, clickedDataId: String(point.z) }),
      )
      setRoisClick(itemId, point.z)
    }
  }

  const onCancel = async () => {
    setAction("")
    if (
      !isEditableRoiPath(refRoiFilePath.current) ||
      workspaceId === undefined
    ) {
      return
    }
    try {
      await dispatch(cancelRoi({ path: refRoiFilePath.current, workspaceId }))
    } finally {
      workspaceId &&
        dispatch(getRoiData({ path: refRoiFilePath.current, workspaceId }))
    }
  }

  useEffect(() => {
    if (!checkStatus && cancelFirst) {
      setCancelFirst(false)
    }
    //eslint-disable-next-line
  }, [JSON.stringify(statusRoi), cancelFirst])

  useEffect(() => {
    if (!roiFilePath || workspaceId === undefined) return
    dispatch(setDataCancel({ roiFilePath, workspaceId, statusRoi }))
    //eslint-disable-next-line
  }, [roiFilePath, action, JSON.stringify(statusRoi)])

  const onAddRoi = () => {
    setAction(ADD_ROI)
  }

  const editRoi = () => {
    setEdit(true)
  }

  const onCancelAdd = () => {
    setAction("")
    setEdit(true)
    setSizeDrag(initSizeDrag)
    setChangeSize(undefined)
  }

  const onMouseDownDragAddRoi = () => {
    setStartDragAddRoi(true)
  }

  const onMouseUpDragAddRoi = () => {
    setStartDragAddRoi(false)
    setChangeSize(undefined)
  }

  const onMouseDownSize = (position: PositionDrag, event: MouseEvent) => {
    setChangeSize(position)
    refPageXSize.current = event.pageX
    refPageYSize.current = event.pageY
  }

  const onMouseMoveAddRoi = (event: MouseEvent<HTMLDivElement>) => {
    const { pageX, pageY } = event
    let newSizeDrag
    if (startDragAddRoi) {
      const { y } = event.currentTarget.getBoundingClientRect()
      let newX = sizeDrag.left + (pageX - refPageXSize.current)
      let newY = Math.ceil(pageY - y - 15) - window.scrollY

      if (newX < 0) newX = 0
      else if (newX + sizeDrag.width > sChart) newX = sChart - sizeDrag.width
      if (newY < 0) newY = 0
      else if (newY + sizeDrag.height > sChart) newY = sChart - sizeDrag.height
      newSizeDrag = { ...sizeDrag, left: newX, top: newY }
    } else if (positionDrag === PositionDrag.LEFT) {
      const newWidth = sizeDrag.width - (pageX - refPageXSize.current)
      const newLeft = sizeDrag.left + (pageX - refPageXSize.current)
      if (newWidth < 10 || newLeft < 1) return
      newSizeDrag = { ...sizeDrag, width: newWidth, left: newLeft }
    } else if (positionDrag === PositionDrag.RIGHT) {
      const newWidth = sizeDrag.width + (pageX - refPageXSize.current)
      if (newWidth < 10 || newWidth > sChart - sizeDrag.left) return
      newSizeDrag = { ...sizeDrag, width: newWidth }
    } else if (positionDrag === PositionDrag.BOTTOM) {
      const newHeight = sizeDrag.height + (pageY - refPageYSize.current)
      if (newHeight < 10 || newHeight > sChart - sizeDrag.top) return
      newSizeDrag = { ...sizeDrag, height: newHeight }
    } else if (positionDrag === PositionDrag.TOP) {
      const newHeight = sizeDrag.height - (pageY - refPageYSize.current)
      const newTop = sizeDrag.top + (pageY - refPageYSize.current)
      if (newHeight < 10 || newTop < 1) return
      newSizeDrag = { ...sizeDrag, height: newHeight, top: newTop }
    }
    if (newSizeDrag) setSizeDrag({ ...sizeDrag, ...newSizeDrag })
    refPageXSize.current = pageX
    refPageYSize.current = pageY
  }

  const addOrSelectRoi = async () => {
    if (!roiFilePath || !workspaceId || !statusRoi) return
    if (action === ADD_ROI) {
      const sizeX = roiDataState[0].length - 1
      const sizeY = roiDataState.length - 1
      const xAdd = Number(((sizeDrag.width + 2) / (sChart / sizeX)).toFixed(1))
      const yAdd = Number(((sizeDrag.height + 2) / (sChart / sizeY)).toFixed(1))
      const x = Number((sizeDrag.left / (sChart / sizeX)).toFixed(1))
      const y = Number((sizeDrag.top / (sChart / sizeY)).toFixed(1))

      const pointCenter = {
        posx: x + Math.floor(xAdd / 2),
        posy: y + Math.floor(yAdd / 2),
        sizex: xAdd,
        sizey: yAdd,
      }
      await dispatch(resetAllOrderList())
      await dispatch(
        addRoi({ path: roiFilePath, workspaceId, data: pointCenter }),
      )
      onCancelAdd()
    }
    if (action === MERGE_ROI) {
      if (roiClicked.length < 2) return
      dispatch(resetAllOrderList())
      await dispatch(
        mergeRoi({
          path: roiFilePath,
          workspaceId,
          data: { ids: roiClicked },
        }),
      )
      resetRoisClick(itemId)
      workspaceId && dispatch(getRoiData({ path: roiFilePath, workspaceId }))
    } else if (action === DELETE_ROI) {
      if (!roiClicked.length) return
      dispatch(resetAllOrderList())
      await dispatch(
        deleteRoi({
          path: roiFilePath,
          workspaceId,
          data: { ids: roiClicked },
        }),
      )
      resetRoisClick(itemId)
      workspaceId && dispatch(getRoiData({ path: roiFilePath, workspaceId }))
    } else if (action === PROMOTE_ROI) {
      if (!roiClicked.length) return
      dispatch(resetAllOrderList())
      try {
        await dispatch(
          promoteRoi({
            path: roiFilePath,
            workspaceId,
            data: { ids: roiClicked },
          }),
        ).unwrap()
      } catch (error) {
        const message = (error as { message?: string })?.message
        enqueueSnackbar(message ?? "Failed to set the ROIs as cell ROIs.", {
          variant: "error",
        })
      }
    }
    setAction("")
    setEdit(true)
    resetRoisClick(itemId)
    dispatch(getStatus({ path: roiFilePath, workspaceId }))
  }

  const onMergeRoi = async () => {
    if (!roiFilePath) return
    setAction(MERGE_ROI)
  }

  const onDeleteRoi = async () => {
    if (!roiFilePath) return
    setAction(DELETE_ROI)
  }

  const onPromoteRoi = async () => {
    if (!roiFilePath) return
    setAction(PROMOTE_ROI)
  }

  const onCommitRoi = async () => {
    if (!roiFilePath || workspaceId === undefined) return
    try {
      await dispatch(commitRoi({ path: roiFilePath, workspaceId })).unwrap()
      // commit rewrites cell_roi and non_cell_roi, so refetch every open view of this node
      const nodeDir = roiFilePath.slice(0, roiFilePath.lastIndexOf("/") + 1)
      const openRoiPaths = Object.values(itemsVisual ?? {})
        .map((item) => (isImageItem(item) ? item.roiItem?.filePath : null))
        .filter(
          (path): path is string =>
            !!path && path !== roiFilePath && path.startsWith(nodeDir),
        )
      await dispatch(getRoiData({ path: roiFilePath, workspaceId })).unwrap()
      openRoiPaths.forEach((path) =>
        dispatch(getRoiData({ path, workspaceId })),
      )

      enqueueSnackbar("Successfully committed to Edit ROI.", {
        variant: "success",
      })
      resetTimeSeries()
      resetRoisClick(itemId)
    } catch (error) {
      enqueueSnackbar("Failed to commit Edit ROI.", { variant: "error" })
    } finally {
      setEdit(false)
      setAction("")
    }
  }

  const resetTimeSeries = () => {
    if (itemsVisual) {
      Object.keys(itemsVisual).forEach((item) => {
        if (isTimeSeriesItem(itemsVisual[item])) {
          dispatch(
            getTimeSeriesInitData({
              path: itemsVisual[item].filePath as string,
              itemId: Number(item),
            }),
          )
        }
      })
    }
  }

  const renderActionRoi = () => {
    if (!allowEditRoi && !allowPromoteRoi) return null
    const okDisabled =
      action !== ADD_ROI && roiClicked.length < (action === MERGE_ROI ? 2 : 1)
    if (action) {
      return (
        <>
          {action !== ADD_ROI ? (
            <BoxDiv data-testid="roi-selected-ids">
              ROI Selecteds: [{roiClicked?.join(",") || ""}]
            </BoxDiv>
          ) : null}
          <BoxDiv sx={{ display: "flex", gap: 1 }}>
            <LinkDiv
              sx={{
                color:
                  action === DELETE_ROI
                    ? "#F84E1B"
                    : action === MERGE_ROI
                      ? "#6619A9"
                      : action === PROMOTE_ROI
                        ? "#32A919"
                        : "default",
                display: "flex",
                gap: 1,
                textDecoration: "none",
                cursor: "default",
              }}
            >
              {`[${action}]`}
            </LinkDiv>
            <LinkDiv onClick={onCancelAdd}>Cancel</LinkDiv>
            <LinkDiv
              style={{
                opacity: okDisabled ? 0.5 : 1,
                cursor: okDisabled ? "default" : "pointer",
              }}
              onClick={addOrSelectRoi}
            >
              OK
            </LinkDiv>
          </BoxDiv>
        </>
      )
    }

    if (!edit)
      return (
        <LinkDiv sx={{ width: "fit-content" }} onClick={editRoi}>
          Edit ROI
        </LinkDiv>
      )
    return null
  }

  if (!statusRoi) return null

  return (
    <ImagePlotContainer>
      <Box sx={{ display: "flex" }}>
        <Box sx={{ flexGrow: 1, mt: 1 }}>
          <PlayBack activeIndex={activeIndex} />
        </Box>
        <FormControlLabel
          sx={{ ml: 1 }}
          control={<Switch checked={selectMode} onChange={handleChange} />}
          label="drag select"
        />
      </Box>
      <Box sx={{ minHeight: 5.5 }}>
        {edit && !action && (allowEditRoi || allowPromoteRoi) ? (
          <>
            <BoxDiv sx={{ flexDirection: "column" }}>
              <BoxWrapper sx={{ marginBottom: 2 }}>
                {allowEditRoi ? (
                  <>
                    <LinkDiv onClick={onAddRoi}>{ADD_ROI}</LinkDiv>
                    <LinkDiv
                      sx={{
                        color: "#F84E1B",
                      }}
                      onClick={onDeleteRoi}
                    >
                      {DELETE_ROI}
                    </LinkDiv>
                    <LinkDiv
                      sx={{
                        color: "#6619A9",
                        ml: 0,
                      }}
                      onClick={onMergeRoi}
                    >
                      {MERGE_ROI}
                    </LinkDiv>
                  </>
                ) : (
                  <LinkDiv
                    data-testid="roi-promote"
                    sx={{
                      color: "#32A919",
                    }}
                    onClick={onPromoteRoi}
                  >
                    {PROMOTE_ROI}
                  </LinkDiv>
                )}
              </BoxWrapper>
              <BoxWrapper>
                {Object.keys(statusRoi).some(
                  (key) => statusRoi[key as keyof StatusROI].length > 0,
                ) ? (
                  <LinkDiv
                    sx={{
                      color: "#32A919",
                    }}
                    onClick={onCommitRoi}
                  >
                    Commit Edit
                  </LinkDiv>
                ) : null}
                <LinkDiv
                  onClick={() => {
                    onCancel()
                    setEdit(false)
                  }}
                >
                  Cancel
                </LinkDiv>
              </BoxWrapper>
            </BoxDiv>
          </>
        ) : (
          renderActionRoi()
        )}
      </Box>
      <div style={{ position: "relative" }}>
        <PlotlyChart
          data={data}
          layout={layout}
          config={config}
          onClick={onChartClick}
          onSelecting={onSelecting}
        />
        {action === ADD_ROI ? (
          <DivAddRoi>
            <DivSvg
              onMouseLeave={onMouseUpDragAddRoi}
              onMouseMove={onMouseMoveAddRoi}
              onMouseUp={onMouseUpDragAddRoi}
            >
              <DivDrag style={sizeDrag}>
                <DragCenter
                  onMouseDown={onMouseDownDragAddRoi}
                  style={{
                    width: sizeDrag.width - 1,
                    height: sizeDrag.height - 1,
                    cursor: !startDragAddRoi ? "grab" : "grabbing",
                  }}
                />
                <DragSizeLeft
                  onMouseDown={(event) =>
                    onMouseDownSize(PositionDrag.LEFT, event)
                  }
                />
                <DragSizeRight
                  onMouseDown={(event) => {
                    onMouseDownSize(PositionDrag.RIGHT, event)
                  }}
                />
                <DragSizeTop
                  onMouseDown={(event) => {
                    onMouseDownSize(PositionDrag.TOP, event)
                  }}
                />
                <DragSizeBottom
                  onMouseDown={(event) => {
                    onMouseDownSize(PositionDrag.BOTTOM, event)
                  }}
                />
              </DivDrag>
            </DivSvg>
          </DivAddRoi>
        ) : null}
      </div>
    </ImagePlotContainer>
  )
})

const PlayBack = memo(function PlayBack({ activeIndex }: ActiveIndexProps) {
  const dispatch = useDispatch()
  const { filePath: path, itemId } = useContext(DisplayDataContext)

  const maxSize = useSelector(selectImageDataMaxSize(path))
  const startIndex = useSelector(selectImageItemStartIndex(itemId))
  const endIndex = useSelector(selectImageItemEndIndex(itemId))
  const duration = useSelector(selectImageItemDuration(itemId))

  const onSliderChange = (event: Event, value: number | number[]) => {
    if (typeof value === "number") {
      const newIndex = value - startIndex
      if (newIndex >= 0 && newIndex !== activeIndex) {
        dispatch(setImageActiveIndex({ itemId, activeIndex: newIndex }))
      }
    }
  }

  const intervalRef = useRef<null | NodeJS.Timeout>(null)

  useEffect(() => {
    if (intervalRef.current !== null) {
      if (activeIndex >= maxSize) {
        clearInterval(intervalRef.current)
        intervalRef.current = null
      }
    }
  }, [activeIndex, maxSize])

  const onPlayClick = useCallback(() => {
    if (activeIndex >= maxSize) {
      dispatch(setImageActiveIndex({ itemId, activeIndex: 0 }))
    }
    if (maxSize > 1 && intervalRef.current === null) {
      intervalRef.current = setInterval(() => {
        dispatch(incrementImageActiveIndex({ itemId }))
      }, duration)
    }
  }, [activeIndex, maxSize, dispatch, duration, itemId])

  const onPauseClick = () => {
    if (intervalRef.current !== null) {
      clearInterval(intervalRef.current)
      intervalRef.current = null
    }
  }

  const onDurationChange = useCallback(
    (event: ChangeEvent<HTMLInputElement>) => {
      const newValue =
        event.target?.value === "" ? "" : Number(event.target?.value)
      if (typeof newValue === "number") {
        dispatch(setImageItemDuration({ itemId, duration: newValue }))
      }
    },
    [dispatch, itemId],
  )
  return (
    <>
      <Button sx={{ mt: 1.5 }} variant="outlined" onClick={onPlayClick}>
        Play
      </Button>
      <Button sx={{ mt: 1.5, ml: 1 }} variant="outlined" onClick={onPauseClick}>
        Pause
      </Button>
      <TextField
        sx={{ width: 100, ml: 2 }}
        label="msec/frame"
        type="number"
        inputProps={{
          step: 100,
          min: 0,
          max: 1000,
        }}
        InputLabelProps={{
          shrink: true,
        }}
        onChange={onDurationChange}
        value={duration}
      />
      <Slider
        aria-label="Custom marks"
        defaultValue={20}
        value={startIndex + activeIndex}
        valueLabelDisplay="auto"
        step={1}
        marks
        min={startIndex}
        max={maxSize === 0 ? 0 : endIndex}
        onChange={onSliderChange}
      />
    </>
  )
})

function imageDataEqualtyFn(
  a: number[][] | undefined,
  b: number[][] | undefined,
) {
  if (a != null && b != null) {
    return twoDimarrayEqualityFn(a, b)
  } else {
    return a === undefined && b === undefined
  }
}

interface PlotDatum {
  curveNumber: number
  data: PlotData
  customdata: Datum
  pointIndex: number
  pointNumber: number
  x: Datum
  xaxis: LayoutAxis
  y: Datum
  yaxis: LayoutAxis
  z: number
}

export function rgba2hex(rgba: number[], alpha: number) {
  const r = rgba[0]
  const g = rgba[1]
  const b = rgba[2]
  const a = alpha

  const outParts = [
    r.toString(16),
    g.toString(16),
    b.toString(16),
    Math.round(a * 255)
      .toString(16)
      .substring(0, 2),
  ]

  // Pad single-digit output values
  outParts.forEach(function (part, i) {
    if (part.length === 1) {
      outParts[i] = "0" + part
    }
  })

  return `#${outParts.join("")}`
}

function debounce<T extends (...args: PlotSelectionEvent[]) => unknown>(
  callback: T,
  delay = 500,
): (...args: Parameters<T>) => void {
  let timeoutId: NodeJS.Timeout
  return (...args) => {
    clearTimeout(timeoutId)
    timeoutId = setTimeout(() => callback(...args), delay)
  }
}

const ImagePlotContainer = styled("div")`
  .selectionlayer path {
    stroke: white !important;
    stroke-dasharray: 6 !important;
    stroke-width: 2px !important;
    opacity: 1 !important;
  }
`

const BoxDiv = styled("div")({
  mt: 1,
  display: "flex",
  listStyle: "none",
  padding: 0,
  margin: 0,
})

const BoxWrapper = styled("div")({
  display: "flex",
  gap: 12,
})

const LinkDiv = styled("div")({
  textDecoration: "underline",
  cursor: "pointer",
  color: "#1155cc",
  zIndex: 999,
  position: "relative",
})

const DivAddRoi = styled("div")({
  width: "100%",
  height: "100%",
  position: "absolute",
  left: 0,
  top: 0,
  borderRadius: 100,
})

const DivSvg = styled("div")({
  width: 321,
  height: 321,
  marginTop: 30,
  marginLeft: 99,
  position: "relative",
})

const DivDrag = styled("div")({
  border: "1px solid #ffffff",
  position: "absolute",
  borderRadius: 100,
})

const DragCenter = styled("div")({
  borderRadius: 100,
  cursor: "grab",
})

const DragSize = styled("div")({
  width: 3,
  height: 3,
  borderRadius: 100,
  position: "absolute",
  background: "#fff",
})

const DragSizeLeft = styled(DragSize)({
  top: "calc(50% - 1px)",
  left: -2,
  cursor: "ew-resize",
})

const DragSizeRight = styled(DragSize)({
  top: "calc(50% - 1px)",
  right: -2,
  cursor: "ew-resize",
})

const DragSizeTop = styled(DragSize)({
  top: -2,
  right: "calc(50% - 1px)",
  cursor: "ns-resize",
})

const DragSizeBottom = styled(DragSize)({
  bottom: -2,
  right: "calc(50% - 1px)",
  cursor: "ns-resize",
})
