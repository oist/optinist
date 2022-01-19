import React from 'react'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

export const HeatmapItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  return <div>HeatmapItemEditor</div>
}
