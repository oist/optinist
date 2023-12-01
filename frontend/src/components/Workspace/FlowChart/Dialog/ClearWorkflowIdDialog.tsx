import { memo } from "react"

import {
  ConfirmDialog,
  ConfirmDialogProps,
} from "components/common/ConfirmDialog"

export const ClearWorkflowIdDialog = memo(function ClearWorkflowIdDialog({
  open,
  onConfirm,
  onCancel,
}: Pick<ConfirmDialogProps, "open" | "onConfirm" | "onCancel">) {
  return (
    <ConfirmDialog
      open={open}
      onConfirm={onConfirm}
      onCancel={onCancel}
      title="Change this parameter?"
      content={`Changing this parameter may affect whole workflow.
      So, this workflow will be run as new workflow with new ID (RUN ALL).
      Existing nodes are kept as they are.`}
      iconType="warning"
    />
  )
})
