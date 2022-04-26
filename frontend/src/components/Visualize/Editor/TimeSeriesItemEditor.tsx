import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import {
  AccordionDetails,
  AccordionSummary,
  FormControlLabel,
  Switch,
  TextField,
} from '@mui/material'
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
  selectTimeSeriesItemFilePath,
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
import {
  getTimeSeriesAllData,
  getTimeSeriesDataById,
} from 'store/slice/DisplayData/DisplayDataActions'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { Accordion } from 'components/Accordion'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'

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
            inputProps={{
              step: 1,
              min: 0,
            }}
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
  const checkedList = useSelector(
    selectTimeSeriesItemCheckedList(itemId),
    // arrayEqualityFn,
  )
  const displayNumbers = useSelector(selectTimeSeriesItemDisplayNumbers(itemId))
  const filePath = useSelector(selectTimeSeriesItemFilePath(itemId))

  const allHandleChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    dispatch(
      setTimeSeriesItemCheckedList({
        itemId,
        checkedList: Object.keys(checkedList).map((_) => {
          return event.target.checked
        }),
      }),
    )

    if (event.target.checked) {
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: Object.keys(checkedList).map((i, v) => {
            return v
          }),
        }),
      )
      if (filePath !== null) {
        dispatch(getTimeSeriesAllData({ path: filePath }))
      }
    } else {
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: [],
        }),
      )
    }
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
        checkedList: Object.values(checkedList).map((v, i) => {
          if (i === index) {
            return event.target.checked
          }
          return v
        }),
      }),
    )

    if (filePath !== null) {
      dispatch(getTimeSeriesDataById({ path: filePath, index }))
    }
  }

  const children = (
    <Box sx={{ display: 'flex', flexDirection: 'column', ml: 3 }}>
      {Object.entries(checkedList).map(([key, value]) => (
        <FormControlLabel
          key={`${key}`}
          label={`Index ${parseInt(key) + 1}`}
          control={
            <Checkbox checked={value} onChange={handleChange} value={key} />
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
                checked={Object.values(checkedList).every((v) => {
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
