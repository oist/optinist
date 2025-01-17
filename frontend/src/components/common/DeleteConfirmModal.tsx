import { FC, useState } from "react"

import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  styled,
  Typography,
} from "@mui/material"

import Input from "components/common/Input"
import Loading from "components/common/Loading"

type DeleteConfirmModalProps = {
  onClose: () => void
  open: boolean
  onSubmit: () => void
  titleSubmit: string
  description: string
  loading?: boolean
}
const DeleteConfirmModal: FC<DeleteConfirmModalProps> = ({
  onClose,
  open,
  onSubmit,
  loading,
  titleSubmit,
  description,
}) => {
  const [textDelete, setTextDelete] = useState("")

  const onConfirm = () => {
    if (textDelete !== "DELETE") return
    onSubmit?.()
    setTextDelete("")
  }

  return (
    <>
      <Dialog open={open} onClose={onClose} maxWidth={"xs"}>
        <DialogTitle>
          <Typography style={{ whiteSpace: "pre-wrap" }}>
            {description}
            This operation cannot be undone. To continue, type
            <span style={{ fontWeight: 600 }}>DELETE</span> in the box below:
          </Typography>
        </DialogTitle>
        <DialogContent>
          <BoxConfirm>
            <Input
              placeholder="DELETE"
              value={textDelete}
              onChange={(e) => setTextDelete(e.target.value)}
              sx={{ width: "calc(100% - 20px)" }}
            />
          </BoxConfirm>
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose} variant={"outlined"}>
            CANCEL
          </Button>
          <Button onClick={onConfirm} color={"error"} variant="contained">
            {titleSubmit}
          </Button>
        </DialogActions>
      </Dialog>
      <Loading loading={loading} />
    </>
  )
}

const BoxConfirm = styled(Box)({
  margin: "20px 0 0",
})

export default DeleteConfirmModal
