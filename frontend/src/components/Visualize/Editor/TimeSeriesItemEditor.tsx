import React from 'react'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

export const TimeSeriesItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  return <div>TimeSeriesItemEditor</div>
}
