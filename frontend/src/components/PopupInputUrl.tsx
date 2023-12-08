import { ChangeEvent, useState } from "react"

import Button from "@mui/material/Button"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogContent from "@mui/material/DialogContent"
import DialogTitle from "@mui/material/DialogTitle"
import TextField from "@mui/material/TextField"
import Typography from "@mui/material/Typography"

import { ConfirmDialog } from "components/common/ConfirmDialog"

type PopupInputUrlProps = {
  open: boolean
  value: string
  setValue: (value: string) => void
  handleClose: () => void
  onLoadFileViaUrl: () => void
  setError: (value: string) => void
  error: string
}

const PopupInputUrl = ({
  open,
  value,
  setValue,
  handleClose,
  onLoadFileViaUrl,
  setError,
  error,
}: PopupInputUrlProps) => {
  const [openConfirm, setOpenConfirm] = useState(false)
  const validateUrl = (url: string) => {
    const validExtensions = [".tiff", ".tif", ".csv", ".hdf5", ".nwb"]
    const fileExtension = url.substring(url.lastIndexOf("."))
    return (
      validExtensions.includes(fileExtension) &&
      (url.startsWith("http://") || url.startsWith("https://"))
    )
  }

  const handleFileVia = (event: ChangeEvent<HTMLInputElement>) => {
    const { value } = event.target
    setValue(value)
    if (!validateUrl(value)) {
      setError("is validate")
    } else {
      setError("")
    }
  }

  const handleClosePopup = () => {
    setError("")
    handleClose()
  }

  const onBlur = () => {
    if (!validateUrl(value)) {
      setError("is validate")
    } else {
      setError("")
    }
  }

  const onClickLoad = () => {
    setOpenConfirm(true)
  }

  return (
    <>
      <Dialog open={open} onClose={handleClosePopup} fullWidth>
        <DialogTitle>Upload file via URL</DialogTitle>
        <DialogContent>
          <TextField
            margin="dense"
            id="name"
            label="Link url"
            fullWidth
            variant="standard"
            value={value}
            onChange={handleFileVia}
            onBlur={onBlur}
          />
          <Typography sx={{ color: "red" }}>{error}</Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleClosePopup} variant={"outlined"}>
            Cancel
          </Button>
          <Button onClick={onClickLoad} variant={"contained"}>
            Load
          </Button>
        </DialogActions>
      </Dialog>
      <ConfirmDialog
        open={openConfirm}
        content={"Do you want load from via url"}
        onCancel={() => setOpenConfirm(false)}
        onConfirm={() => {
          setOpenConfirm(false)
          onLoadFileViaUrl()
          setValue("")
        }}
      />
    </>
  )
}

export default PopupInputUrl
