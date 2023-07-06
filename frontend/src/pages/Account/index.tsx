import { useSelector, useDispatch } from 'react-redux'
import { Box, styled, Typography } from '@mui/material'
import Loading from "components/common/Loading"
import ChangePasswordModal from 'components/Account/ChangePasswordModal'
import DeleteConfirmModal from 'components/common/DeleteConfirmModal'
import { useState } from 'react'
import { useNavigate } from "react-router-dom";
import { updateMePasswordApi } from 'api/users/UsersMe'
import { deleteMe } from 'store/slice/User/UserActions'
import { selectCurrentUser } from 'store/slice/User/UserSelector'
const Account = () => {
  const user = useSelector(selectCurrentUser)
  const dispatch = useDispatch()
  const navigate = useNavigate()
  const [isDeleteConfirmModalOpen, setIsDeleteConfirmModalOpen] = useState(false)
  const [isChangePwModalOpen, setIsChangePwModalOpen] = useState(false)
  const [isLoading, setIsLoading] = useState(false)

  const handleCloseDeleteComfirmModal = () => {
    setIsDeleteConfirmModalOpen(false)
  }

  const onDeleteAccountClick = () => {
    setIsDeleteConfirmModalOpen(true)
  }

  const onConfirmDelete = async () => {
    if(!user) return
    setIsLoading(true)
    try {
      dispatch(deleteMe())
      navigate('/login')
    }
    catch {}
    finally {
      setIsLoading(false)
    }
    handleCloseDeleteComfirmModal()
  }

  const handleCloseChangePw = () => {
    setIsChangePwModalOpen(false)
  }

  const onChangePwClick = () => {
    setIsChangePwModalOpen(true)
  }

  const onConfirmChangePw = async (oldPass: string, newPass: string) => {
    setIsLoading(true)
    try {
      await updateMePasswordApi({old_password: oldPass, new_password: newPass})
      alert('Your password has been successfully changed.')
      handleCloseChangePw()
    }
    catch {
      alert('Failed to Change Password!')
    }
    finally {
      setIsLoading(false)
    }
  }

  return (
    <AccountWrapper>
      <DeleteConfirmModal
        titleSubmit="Delete My Account"
        onClose={handleCloseDeleteComfirmModal}
        open={isDeleteConfirmModalOpen}
        onSubmit={onConfirmDelete}
        description='Delete account will erase all of your data.'
      />
      <ChangePasswordModal
        onSubmit={onConfirmChangePw}
        open={isChangePwModalOpen}
        onClose={handleCloseChangePw}
      />
      <Title>Account Profile</Title>
      <BoxFlex>
        <TitleData>Account ID</TitleData>
        <BoxData>{user?.uid}</BoxData>
      </BoxFlex>
      <BoxFlex>
        <TitleData>Email</TitleData>
        <BoxData>{user?.email}</BoxData>
      </BoxFlex>
      <BoxFlex sx={{ justifyContent: 'space-between', mt: 10 }}>
        <ButtonSubmit onClick={onChangePwClick}>Change Password</ButtonSubmit>
        <ButtonSubmit onClick={onDeleteAccountClick}>Delete Account</ButtonSubmit>
      </BoxFlex>
      {
        isLoading && <Loading />
      }
    </AccountWrapper>
  )
}

const AccountWrapper = styled(Box)({
  padding: '0 20px',
})

const BoxFlex = styled(Box)({
  display: 'flex',
  margin: '20px 0 10px 0',
  alignItems: 'center',
  maxWidth: 1000,
})

const Title = styled('h2')({
  marginBottom: 40,
})

const BoxData = styled(Typography)({
  fontWeight: 700,
  minWidth: 272,
})

const TitleData = styled(Typography)({
  width: 250,
})

const ButtonSubmit = styled('button')({
  backgroundColor: '#283237',
  color: '#ffffff',
  borderRadius: 4,
  border: 'none',
  outline: 'none',
  padding: '10px 20px',
  cursor: 'pointer',
})

export default Account
