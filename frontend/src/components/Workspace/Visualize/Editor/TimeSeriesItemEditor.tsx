import { ChangeEvent, FC, useContext, useEffect, useState } from "react"
import { useSelector, useDispatch } from "react-redux"

import ExpandMoreIcon from "@mui/icons-material/ExpandMore"
import {
  AccordionDetails,
  AccordionSummary,
  FormControlLabel,
  Grid,
  MenuItem,
  SelectChangeEvent,
} from "@mui/material"
import Box from "@mui/material/Box"
import Checkbox from "@mui/material/Checkbox"

import { Accordion } from "components/common/Accordion"
import { ParamSection } from "components/common/ParamSection"
import { ParamSelect } from "components/common/ParamSelect"
import { ParamSwitch } from "components/common/ParamSwitch"
import { ParamTextField } from "components/common/ParamTextField"
import { SaveFig } from "components/Workspace/Visualize/Editor/SaveFig"
import { SelectedItemIdContext } from "components/Workspace/Visualize/VisualizeItemEditor"
import {
  getTimeSeriesAllData,
  getTimeSeriesDataById,
} from "store/slice/DisplayData/DisplayDataActions"
import { selectFrameRate } from "store/slice/Experiments/ExperimentsSelectors"
import { selectPipelineLatestUid } from "store/slice/Pipeline/PipelineSelectors"
import {
  selectTimeSeriesItemDrawOrderList,
  selectTimeSeriesItemOffset,
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
  selectTimeSeriesItemFilePath,
  selectTimeSeriesItemKeys,
  selectImageItemRangeUnit,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"
import {
  setTimeSeriesItemOffset,
  setTimeSeriesItemShowGrid,
  setTimeSeriesItemShowLine,
  setTimeSeriesItemShowTickLabels,
  setTimeSeriesItemSpan,
  setTimeSeriesItemXrangeLeft,
  setTimeSeriesItemXrangeRight,
  setTimeSeriesItemZeroLine,
  setTimeSeriesItemDrawOrderList,
  changeRangeUnit,
  setClickedData,
} from "store/slice/VisualizeItem/VisualizeItemSlice"
import { AppDispatch } from "store/store"
import { arrayEqualityFn } from "utils/EqualityUtils"

export const TimeSeriesItemEditor: FC = () => {
  return (
    <>
      <ParamSection title="TimeSeries">
        <STD />
        <Span />
        <ShowGrid />
        <ShowLine />
        <ShowTickLabels />
        <ZeroLine />
        <SelectValue />
        <Xrange />
        <LegendSelect />
      </ParamSection>
      <SaveFig />
    </>
  )
}

const STD: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const stdBool = useSelector(selectTimeSeriesItemOffset(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemOffset({ itemId, stdBool: !stdBool }))
  }
  return <ParamSwitch label="STD" value={stdBool} onChange={toggleChecked} />
}

const Span: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const span = useSelector(selectTimeSeriesItemSpan(itemId))

  const dispatch = useDispatch()
  const onChange = (event: ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newValue === "number" && newValue > 0) {
      dispatch(setTimeSeriesItemSpan({ itemId, span: newValue }))
    }
  }
  return (
    <ParamTextField
      type="number"
      label="vertical offset"
      value={span ?? ""}
      onChange={onChange}
    />
  )
}

const ShowGrid: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showgrid = useSelector(selectTimeSeriesItemShowGrid(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemShowGrid({ itemId, showgrid: !showgrid }))
  }
  return (
    <ParamSwitch label="ShowGrid" value={showgrid} onChange={toggleChecked} />
  )
}

const ShowLine: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showline = useSelector(selectTimeSeriesItemShowLine(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemShowLine({ itemId, showline: !showline }))
  }
  return (
    <ParamSwitch label="ShowLine" value={showline} onChange={toggleChecked} />
  )
}

const ShowTickLabels: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const showticklabels = useSelector(selectTimeSeriesItemShowTickLabels(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(
      setTimeSeriesItemShowTickLabels({
        itemId,
        showticklabels: !showticklabels,
      }),
    )
  }
  return (
    <ParamSwitch
      label="ShowTickLabels"
      value={showticklabels}
      onChange={toggleChecked}
    />
  )
}

const ZeroLine: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const zeroline = useSelector(selectTimeSeriesItemZeroLine(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemZeroLine({ itemId, zeroline: !zeroline }))
  }
  return (
    <ParamSwitch label={"ZeroLine"} value={zeroline} onChange={toggleChecked} />
  )
}

const SelectValue: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const value = useSelector(selectImageItemRangeUnit(itemId))
  const dispatch = useDispatch()
  const onChangeValue = async (e: SelectChangeEvent) => {
    dispatch(changeRangeUnit({ itemId: itemId, rangeUnit: e.target.value }))
  }
  return (
    <ParamSelect label={"Range Unit"} value={value} onChange={onChangeValue}>
      <MenuItem value={"frames"}>Frames</MenuItem>
      <MenuItem value={"time"}>Time</MenuItem>
    </ParamSelect>
  )
}

const Xrange: FC = () => {
  const currentPipelineUid = useSelector(selectPipelineLatestUid)
  const itemId = useContext(SelectedItemIdContext)
  const frameRate = useSelector(selectFrameRate(currentPipelineUid))
  const rangeUnit = useSelector(selectImageItemRangeUnit(itemId))
  const xrangeSelector = useSelector(selectTimeSeriesItemXrange(itemId))

  const [xrange, setXrange] = useState(xrangeSelector)

  useEffect(() => {
    if (Object.keys(xrange).length < 1) return
    rangeUnit === "frames"
      ? setXrange(xrangeSelector)
      : setXrange({
          left: Number(xrangeSelector.left) / frameRate,
          right: Number(xrangeSelector.right) / frameRate,
        })
    //eslint-disable-next-line
  }, [JSON.stringify(rangeUnit), JSON.stringify(xrangeSelector)])

  const dispatch = useDispatch()
  const onChangeLeft = (event: ChangeEvent<HTMLInputElement>) => {
    const newLeft = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newLeft === "number") {
      dispatch(
        setTimeSeriesItemXrangeLeft({
          itemId,
          left: rangeUnit === "frames" ? newLeft : newLeft * frameRate,
        }),
      )
    }
  }
  const onChangeRight = (event: ChangeEvent<HTMLInputElement>) => {
    const newRight = event.target.value === "" ? "" : Number(event.target.value)
    if (typeof newRight === "number") {
      dispatch(
        setTimeSeriesItemXrangeRight({
          itemId,
          right: rangeUnit === "frames" ? newRight : newRight * frameRate,
        }),
      )
    }
  }

  return (
    <Grid container justifyContent="space-between">
      <Grid item>
        <ParamTextField
          label="Left"
          type="number"
          inputProps={{
            min: 0,
          }}
          style={{ width: 105 }}
          onChange={onChangeLeft}
          value={xrange.left ?? ""}
        />
      </Grid>
      <Grid item>
        <ParamTextField
          label="Right"
          type="number"
          style={{ width: 105 }}
          onChange={onChangeRight}
          value={xrange.right ?? ""}
        />
      </Grid>
    </Grid>
  )
}

const LegendSelect: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const dispatch = useDispatch<AppDispatch>()
  // const drawIndexMap = useSelector(selectTimeSeriesItemDrawIndexMap(itemId))
  const dataKeys = useSelector(
    selectTimeSeriesItemKeys(itemId),
    arrayEqualityFn,
  )
  const drawOrderList = useSelector(
    selectTimeSeriesItemDrawOrderList(itemId),
    arrayEqualityFn,
  )
  const filePath = useSelector(selectTimeSeriesItemFilePath(itemId))

  const allHandleChange = (event: ChangeEvent<HTMLInputElement>) => {
    dispatch(
      setTimeSeriesItemDrawOrderList({
        itemId,
        drawOrderList: event.target.checked ? dataKeys : [],
      }),
    )

    if (event.target.checked && filePath !== null) {
      dispatch(getTimeSeriesAllData({ path: filePath }))
    }
  }

  const handleChange = (event: ChangeEvent<HTMLInputElement>) => {
    const index = event.target.value
    const newDrawOrderList = event.target.checked
      ? [...drawOrderList, index]
      : drawOrderList.filter((value) => value !== index)

    dispatch(
      setTimeSeriesItemDrawOrderList({
        itemId,
        drawOrderList: newDrawOrderList,
      }),
    )
    dispatch(
      setClickedData({
        itemId,
        clickedDataId: event.target.checked ? index : null,
      }),
    ),
      // eslint-disable-next-line no-console
      console.log("#############\n In TimeSeriesItemEditor.tsx, clickedDataId")

    if (filePath !== null) {
      dispatch(getTimeSeriesDataById({ path: filePath, index: index }))
    }
  }

  const drawIndexMap = Object.fromEntries(
    dataKeys.map((v) => {
      if (drawOrderList.includes(v)) {
        return [v, true]
      } else {
        return [v, false]
      }
    }),
  )

  const children = (
    <Box sx={{ display: "flex", flexDirection: "column", ml: 3 }}>
      {dataKeys.map((key) => (
        <FormControlLabel
          key={`${key}`}
          label={`Index ${key}`}
          control={
            <Checkbox
              checked={drawIndexMap[key]}
              onChange={handleChange}
              value={key}
            />
          }
        />
      ))}
    </Box>
  )

  return (
    <Accordion sx={{ my: 2 }} TransitionProps={{ unmountOnExit: true }}>
      <AccordionSummary expandIcon={<ExpandMoreIcon />}>
        Legend select
      </AccordionSummary>
      <AccordionDetails>
        <>
          <FormControlLabel
            label="All Check"
            control={
              <Checkbox
                checked={Object.values(drawIndexMap).every((v) => {
                  return v
                })}
                onChange={allHandleChange}
              />
            }
          />
          {children}
        </>
      </AccordionDetails>
    </Accordion>
  )
}
