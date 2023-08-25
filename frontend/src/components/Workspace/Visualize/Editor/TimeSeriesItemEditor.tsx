import React, {useEffect, useState} from 'react'
import {useSelector, useDispatch, DefaultRootState} from 'react-redux'
import {
  AccordionDetails,
  AccordionSummary,
  FormControlLabel, MenuItem, Select,
  Switch,
  TextField,
} from '@mui/material'
import Box from '@mui/material/Box'
import Checkbox from '@mui/material/Checkbox'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'

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
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import {
  setTimeSeriesItemOffset,
  setTimeSeriesItemShowGrid,
  setTimeSeriesItemShowLine,
  setTimeSeriesItemShowTickLabels,
  setTimeSeriesItemSpan,
  setTimeSeriesItemXrangeLeft,
  setTimeSeriesItemXrangeRight,
  setTimeSeriesItemZeroLine,
  setTimeSeriesItemDrawOrderList, changeRangeUnit,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  getTimeSeriesAllData,
  getTimeSeriesDataById,
} from 'store/slice/DisplayData/DisplayDataActions'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { Accordion } from 'components/common/Accordion'

import { SaveFig } from './SaveFig'

export const TimeSeriesItemEditor: React.FC = () => {
  return (
    <div style={{ margin: '10px', padding: 10 }}>
      <Offset />
      <Span />
      <ShowGrid />
      <ShowLine />
      <ShowTickLabels />
      <ZeroLine />
      <SelectValue />
      <Xrange />
      <LegendSelect />
      <SaveFig />
    </div>
  )
}

const Offset: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const offset = useSelector(selectTimeSeriesItemOffset(itemId))
  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemOffset({ itemId, offset: !offset }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={offset} onChange={toggleChecked} />}
      label="offset"
    />
  )
}

const Span: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const span = useSelector(selectTimeSeriesItemSpan(itemId))

  const dispatch = useDispatch()
  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number' && newValue > 0) {
      dispatch(setTimeSeriesItemSpan({ itemId, span: newValue }))
    }
  }
  return (
    <FormControlLabel
      control={
        <TextField
          type="number"
          style={{ width: '6vw' }}
          InputLabelProps={{
            shrink: true,
          }}
          onChange={onChange}
          value={span}
        />
      }
      label="offset std"
    />
  )
}

const ShowGrid: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showgrid = useSelector(selectTimeSeriesItemShowGrid(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemShowGrid({ itemId, showgrid: !showgrid }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={showgrid} onChange={toggleChecked} />}
      label="showgrid"
    />
  )
}

const ShowLine: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const showline = useSelector(selectTimeSeriesItemShowLine(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemShowLine({ itemId, showline: !showline }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={showline} onChange={toggleChecked} />}
      label="showline"
    />
  )
}

const ShowTickLabels: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
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
    <FormControlLabel
      control={<Switch checked={showticklabels} onChange={toggleChecked} />}
      label="showticklabels"
    />
  )
}

const ZeroLine: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const zeroline = useSelector(selectTimeSeriesItemZeroLine(itemId))

  const dispatch = useDispatch()
  const toggleChecked = () => {
    dispatch(setTimeSeriesItemZeroLine({ itemId, zeroline: !zeroline }))
  }
  return (
    <FormControlLabel
      control={<Switch checked={zeroline} onChange={toggleChecked} />}
      label="zeroline"
    />
  )
}

const SelectValue: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const value = useSelector(selectImageItemRangeUnit(itemId))
  const dispatch = useDispatch()
  const onChangeValue = async (e: any) => {
    dispatch(changeRangeUnit({itemId: itemId, rangeUnit: e.target.value}))
  }
  return (
    <Box sx={{marginBottom: 2}}>
      <p>range unit</p>
      <Select
        sx={{width: '100%'}}
        value={value}
        onChange={onChangeValue}
      >
        <MenuItem value={'frames'}>Frames</MenuItem>
        <MenuItem value={'time'}>Time</MenuItem>
      </Select>
    </Box>
  )
}

const Xrange: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)

  const rangeUnit = useSelector(selectImageItemRangeUnit(itemId))
  const xrangeSelector = useSelector(selectTimeSeriesItemXrange(itemId))

  const [xrange, setXrange] = useState(xrangeSelector)

  useEffect(() => {
    if(Object.keys(xrange).length < 1) return
    rangeUnit === 'frames' ? setXrange(xrangeSelector) : setXrange({left: Number(xrangeSelector.left)/ 50, right: Number(xrangeSelector.right) / 50})
  }, [JSON.stringify(rangeUnit), JSON.stringify(xrangeSelector)])

  const dispatch = useDispatch()
  const onChangeLeft = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newLeft = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newLeft === 'number') {
      dispatch(setTimeSeriesItemXrangeLeft({ itemId, left: rangeUnit === 'frames' ? newLeft : newLeft * 50 }))
    }
  }
  const onChangeRight = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newRight = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newRight === 'number') {
      dispatch(setTimeSeriesItemXrangeRight({ itemId, right: rangeUnit === 'frames' ? newRight : newRight * 50 }))
    }
  }

  return (
    <FormControlLabel
      control={
        <>
          <TextField
            label="left"
            style={{ width: 50 }}
            type="number"
            inputProps={{
              step: 1,
              min: 0,
            }}
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onChangeLeft}
            value={xrange.left ?? ''}
            defaultValue={undefined}
          />
          <TextField
            label="right"
            style={{ width: 50 }}
            type="number"
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onChangeRight}
            value={xrange.right ?? ''}
          />
        </>
      }
      label=""
    />
  )
}

const LegendSelect: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const dispatch = useDispatch()
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

  const allHandleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
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

  const handleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
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
    <Box sx={{ display: 'flex', flexDirection: 'column', ml: 3 }}>
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
    <Accordion sx={{ mt: 2 }} TransitionProps={{ unmountOnExit: true }}>
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
