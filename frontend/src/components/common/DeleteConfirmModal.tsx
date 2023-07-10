import { FC, useState } from 'react'
import { Box, Button, Modal, styled, Typography } from '@mui/material'
import Input  from 'components/common/Input'
import Loading from "../common/Loading";

type DeleteConfirmModalProps = {
  onClose: () => void
  open: boolean
  onSubmit: () => void
  titleSubmit: string
  description: string
  isLoading?: boolean
}
const DeleteConfirmModal: FC<DeleteConfirmModalProps> = ({
  onClose,
  open,
  onSubmit,
  isLoading,
  titleSubmit,
  description,
}) => {
  const [textDelete, setTextDelete] = useState('')

  const onConfirm = () => {
    if (textDelete !== 'DELETE') return
    onSubmit?.()
    setTextDelete('')
  }

  return (
    <>
      <Modal
          open={open}
          onClose={onClose}
          aria-labelledby="modal-modal-title"
          aria-describedby="modal-modal-description"
      >
        <ContentDelete>
          <Typography style={{ whiteSpace: 'pre-wrap' }}>
            {description}
            This operation cannot be undone.
            To continue, type "<span style={{ fontWeight: 600 }}>DELETE</span>" in the box below:
          </Typography>
          <BoxConfirm>
            <Input
                placeholder="DELETE"
                value={textDelete}
                onChange={(e) => setTextDelete(e.target.value)}
            />
            <ButtonConfirm onClick={onConfirm} sx={{ backgroundColor: 'red !important' }}>{titleSubmit}</ButtonConfirm>
          </BoxConfirm>
          <Button onClick={onClose}>
            <Typography
                sx={{
                  textDecoration: 'underline',
                  textTransform: 'none',
                  lineHeight: '17px',
                }}
            >
              Close
            </Typography>
          </Button>
        </ContentDelete>
      </Modal>
      {
        isLoading ? <Loading /> : null
      }
    </>
  )
}

const ContentDelete = styled(Box)`
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  width: 500px;
  background-color: rgb(255, 255, 255);
  box-shadow: rgb(0 0 0 / 20%) 0px 11px 15px -7px,
    rgb(0 0 0 / 14%) 0px 24px 38px 3px, rgb(0 0 0 / 12%) 0px 9px 46px 8px;
  padding: 16px;
  border-radius: 4px;
  outline: none;
`

const ButtonConfirm = styled(Button)({
  backgroundColor: '#283237 !important',
  height: 36,
  marginLeft: 10,
  color: '#ffffff',
  marginTop: -1,
})

const BoxConfirm = styled(Box)({
  margin: '20px 0 0',
})

export default DeleteConfirmModal
