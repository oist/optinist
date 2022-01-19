import React from 'react'
import { SelectedItemIdContext } from '../VisualizeItemEditor'

export const TableItemEditor: React.FC = () => {
  const itemId = React.useContext(SelectedItemIdContext)
  return <div>TableItemEditor</div>
}
