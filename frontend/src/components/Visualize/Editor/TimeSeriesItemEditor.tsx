import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { FormControlLabel, Switch, TextField } from '@mui/material'
import Box from '@mui/material/Box'
import Checkbox from '@mui/material/Checkbox'
import {
  selectTimeSeriesItemCheckedList,
  selectTimeSeriesItemDisplayNumbers,
  selectTimeSeriesItemOffset,
  selectTimeSeriesItemShowGrid,
  selectTimeSeriesItemShowLine,
  selectTimeSeriesItemShowTickLabels,
  selectTimeSeriesItemSpan,
  selectTimeSeriesItemXrange,
  selectTimeSeriesItemZeroLine,
  selectVisualizeDataFilePath,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { SelectedItemIdContext } from '../VisualizeItemEditor'
import {
  setTimeSeriesItemCheckedList,
  setTimeSeriesItemDisplayNumbers,
  setTimeSeriesItemOffset,
  setTimeSeriesItemShowGrid,
  setTimeSeriesItemShowLine,
  setTimeSeriesItemShowTickLabels,
  setTimeSeriesItemSpan,
  setTimeSeriesItemXrangeLeft,
  setTimeSeriesItemXrangeRight,
  setTimeSeriesItemZeroLine,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { getTimeSeriesData } from 'store/slice/DisplayData/DisplayDataActions'

export const TimeSeriesItemEditor: React.FC = () => {
  return (
    <div style={{ margin: '10px', padding: 10 }}>
      <Offset />
      <Span />
      <ShowGrid />
      <ShowLine />
      <ShowTickLabels />
      <ZeroLine />
      <Xrange />
      <LegendSelect />
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
    if (typeof newValue === 'number') {
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
          defaultValue={span}
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

const Xrange: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  const xrange = useSelector(selectTimeSeriesItemXrange(itemId))

  const dispatch = useDispatch()
  const onChangeLeft = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newLeft = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newLeft === 'number') {
      dispatch(setTimeSeriesItemXrangeLeft({ itemId, left: newLeft }))
    }
  }
  const onChangeRight = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newRight = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newRight === 'number') {
      dispatch(setTimeSeriesItemXrangeRight({ itemId, right: newRight }))
    }
  }

  return (
    <FormControlLabel
      control={
        <>
          left:
          <TextField
            style={{ width: 50 }}
            type="number"
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onChangeLeft}
            defaultValue={xrange.left}
          />
          right:
          <TextField
            style={{ width: 50 }}
            type="number"
            InputLabelProps={{
              shrink: true,
            }}
            onChange={onChangeRight}
            defaultValue={xrange.right}
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
  const checkedList = useSelector(selectTimeSeriesItemCheckedList(itemId))
  const displayNumbers = useSelector(selectTimeSeriesItemDisplayNumbers(itemId))
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))

  const allHandleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    dispatch(
      setTimeSeriesItemCheckedList({
        itemId,
        checkedList: checkedList.map((_) => {
          return event.target.checked
        }),
      }),
    )
  }

  const handleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const index = parseInt(event.target.value)

    // displayNumbers
    if (event.target.checked) {
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: [...displayNumbers, index],
        }),
      )
    } else {
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: displayNumbers.filter((value) => value !== index),
        }),
      )
    }

    // CheckList
    dispatch(
      setTimeSeriesItemCheckedList({
        itemId,
        checkedList: checkedList.map((v, i) => {
          if (i === index) {
            return event.target.checked
          }
          return v
        }),
      }),
    )

    if (filePath !== null) {
      dispatch(getTimeSeriesData({ path: filePath, index }))
    }
  }

  const children = (
    <Box sx={{ display: 'flex', flexDirection: 'column', ml: 3 }}>
      {checkedList.map((v, i) => (
        <FormControlLabel
          key={`${i}`}
          label={`Index ${i + 1}`}
          control={<Checkbox checked={v} onChange={handleChange} value={i} />}
        />
      ))}
    </Box>
  )

  return (
    <div>
      <FormControlLabel
        label="All Check"
        control={
          <Checkbox
            checked={checkedList.every((v) => {
              return v
            })}
            onChange={allHandleChange}
          />
        }
      />
      {children}
    </div>
  )
}
