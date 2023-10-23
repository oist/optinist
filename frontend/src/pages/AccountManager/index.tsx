import {
  ChangeEvent,
  useCallback,
  useEffect,
  useMemo,
  useState,
  MouseEvent,
} from "react"
import { useDispatch, useSelector } from "react-redux"
import { useNavigate, useSearchParams } from "react-router-dom"

import { useSnackbar, VariantType } from "notistack"

import DeleteIcon from "@mui/icons-material/Delete"
import EditIcon from "@mui/icons-material/Edit"
import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogTitle,
  Input,
  styled,
} from "@mui/material"
import IconButton from "@mui/material/IconButton"
import { SelectChangeEvent } from "@mui/material/Select"
import {
  DataGrid,
  GridColDef,
  GridFilterInputValueProps,
  GridFilterModel,
  GridRenderCellParams,
  GridSortDirection,
  GridSortItem,
  GridSortModel,
} from "@mui/x-data-grid"
import { isRejectedWithValue } from "@reduxjs/toolkit"

import { ROLE } from "@types"
import { AddUserDTO, UserDTO } from "api/users/UsersApiDTO"
import InputError from "components/common/InputError"
import Loading from "components/common/Loading"
import PaginationCustom from "components/common/PaginationCustom"
import SelectError from "components/common/SelectError"
import { regexEmail, regexIgnoreS, regexPassword } from "const/Auth"
import {
  deleteUser,
  createUser,
  getListUser,
  updateUser,
} from "store/slice/User/UserActions"
import {
  isAdmin,
  selectCurrentUser,
  selectListUser,
  selectLoading,
} from "store/slice/User/UserSelector"
import { AppDispatch } from "store/store"

let timeout: NodeJS.Timeout | undefined = undefined

type ModalComponentProps = {
  onSubmitEdit: (
    id: number | string | undefined,
    data: { [key: string]: string },
  ) => void
  setOpenModal: (v: boolean) => void
  dataEdit?: {
    [key: string]: string
  }
}

type PopupType = {
  open: boolean
  handleClose: () => void
  handleOkDel: () => void
  name?: string
}

const initState = {
  email: "",
  password: "",
  role_id: "",
  name: "",
  confirmPassword: "",
}

const ModalComponent = ({
  onSubmitEdit,
  setOpenModal,
  dataEdit,
}: ModalComponentProps) => {
  const [formData, setFormData] = useState<{ [key: string]: string }>(
    dataEdit || initState,
  )
  const [isDisabled, setIsDisabled] = useState(false)
  const [errors, setErrors] = useState<{ [key: string]: string }>(initState)

  const validateEmail = (value: string): string => {
    const error = validateField("email", 255, value)
    if (error) return error
    if (!regexEmail.test(value)) {
      return "Invalid email format"
    }
    return ""
  }

  const validatePassword = (
    value: string,
    isConfirm = false,
    values?: { [key: string]: string },
  ): string => {
    if (!value && !dataEdit?.uid) return "This field is required"
    const errorLength = validateLength("password", 255, value)
    if (errorLength) {
      return errorLength
    }
    const datas = values || formData
    if (!regexPassword.test(value) && value) {
      return "Your password must be at least 6 characters long and must contain at least one letter, number, and special character"
    }
    if (regexIgnoreS.test(value)) {
      return "Allowed special characters (!#$%&()*+,-./@_|)"
    }
    if (isConfirm && datas.password !== value && value) {
      return "password is not match"
    }
    return ""
  }

  const validateField = (name: string, length: number, value?: string) => {
    if (!value) return "This field is required"
    return validateLength(name, length, value)
  }

  const validateLength = (name: string, length: number, value?: string) => {
    if (value && value.length > length)
      return `The text may not be longer than ${length} characters`
    if (formData[name]?.length && value && value.length > length) {
      return `The text may not be longer than ${length} characters`
    }
    return ""
  }

  const validateForm = (): { [key: string]: string } => {
    const errorName = validateField("name", 100, formData.name)
    const errorEmail = validateEmail(formData.email)
    const errorRole = validateField("role_id", 50, formData.role_id)
    const errorPassword = dataEdit?.id
      ? ""
      : validatePassword(formData.password)
    const errorConfirmPassword = dataEdit?.id
      ? ""
      : validatePassword(formData.confirmPassword, true)
    return {
      email: errorEmail,
      password: errorPassword,
      confirmPassword: errorConfirmPassword,
      name: errorName,
      role_id: errorRole,
    }
  }

  const onChangeData = (
    e: ChangeEvent<HTMLTextAreaElement | HTMLInputElement> | SelectChangeEvent,
    length: number,
  ) => {
    const { value, name } = e.target
    const newDatas = { ...formData, [name]: value }
    setFormData(newDatas)
    let error: string =
      name === "email"
        ? validateEmail(value)
        : validateField(name, length, value)
    let errorConfirm = errors.confirmPassword
    if (name.toLowerCase().includes("password")) {
      error = validatePassword(value, name === "confirmPassword", newDatas)
      if (name !== "confirmPassword" && formData.confirmPassword) {
        errorConfirm = validatePassword(
          newDatas.confirmPassword,
          true,
          newDatas,
        )
      }
    }
    setErrors({ ...errors, confirmPassword: errorConfirm, [name]: error })
  }

  const onSubmit = async (e: MouseEvent<HTMLButtonElement>) => {
    e.preventDefault()
    setIsDisabled(true)
    const newErrors = validateForm()
    if (Object.keys(newErrors).some((key) => !!newErrors[key])) {
      setErrors(newErrors)
      setIsDisabled(false)
      return
    }
    try {
      await onSubmitEdit(dataEdit?.id, formData)
      setOpenModal(false)
    } finally {
      setIsDisabled(false)
    }
  }
  const onCancel = () => {
    setOpenModal(false)
  }

  return (
    <Modal>
      <ModalBox>
        <TitleModal>{dataEdit?.id ? "Edit" : "Add"} Account</TitleModal>
        <BoxData>
          <LabelModal>Name: </LabelModal>
          <InputError
            name="name"
            value={formData?.name || ""}
            onChange={(e) => onChangeData(e, 100)}
            onBlur={(e) => onChangeData(e, 100)}
            errorMessage={errors.name}
          />
          <LabelModal>Role: </LabelModal>
          <SelectError
            value={formData?.role_id || ""}
            options={Object.keys(ROLE).filter((key) => !Number(key))}
            name="role_id"
            onChange={(e) => onChangeData(e, 50)}
            onBlur={(e) => onChangeData(e, 50)}
            errorMessage={errors.role_id}
          />
          <LabelModal>e-mail: </LabelModal>
          <InputError
            name="email"
            value={formData?.email || ""}
            onChange={(e) => onChangeData(e, 255)}
            onBlur={(e) => onChangeData(e, 255)}
            errorMessage={errors.email}
          />
          {!dataEdit?.id ? (
            <>
              <LabelModal>Password: </LabelModal>
              <InputError
                name="password"
                value={formData?.password || ""}
                onChange={(e) => onChangeData(e, 255)}
                onBlur={(e) => onChangeData(e, 255)}
                type={"password"}
                errorMessage={errors.password}
              />
              <LabelModal>Confirm Password: </LabelModal>
              <InputError
                name="confirmPassword"
                value={formData?.confirmPassword || ""}
                onChange={(e) => onChangeData(e, 255)}
                onBlur={(e) => onChangeData(e, 255)}
                type={"password"}
                errorMessage={errors.confirmPassword}
              />
            </>
          ) : null}
        </BoxData>
        <ButtonModal>
          <Button onClick={() => onCancel()}>Cancel</Button>
          <Button disabled={isDisabled} onClick={(e) => onSubmit(e)}>
            Ok
          </Button>
        </ButtonModal>
      </ModalBox>
      {isDisabled ? <Loading /> : null}
    </Modal>
  )
}

const PopupDelete = ({ open, handleClose, handleOkDel, name }: PopupType) => {
  if (!open) return null
  return (
    <Box>
      <Dialog open={open} onClose={handleClose} sx={{ margin: 0 }}>
        <DialogTitle>{`Do you want delete User "${name}"?`}</DialogTitle>
        <DialogActions>
          <Button onClick={handleClose}>Cancel</Button>
          <Button onClick={handleOkDel}>Ok</Button>
        </DialogActions>
      </Dialog>
    </Box>
  )
}

const AccountManager = () => {
  const dispatch = useDispatch<AppDispatch>()
  const navigate = useNavigate()

  const listUser = useSelector(selectListUser)
  const loading = useSelector(selectLoading)
  const user = useSelector(selectCurrentUser)
  const admin = useSelector(isAdmin)

  const [openModal, setOpenModal] = useState(false)
  const [dataEdit, setDataEdit] = useState({})
  const [newParams, setNewParams] = useState(
    window.location.search.replace("?", ""),
  )
  const [openDel, setOpenDel] = useState<{
    id?: number
    name?: string
    open: boolean
  }>()

  const [searchParams, setParams] = useSearchParams()
  const limit = searchParams.get("limit") || 50
  const offset = searchParams.get("offset") || 0
  const name = searchParams.get("name") || undefined
  const email = searchParams.get("email") || undefined
  const sort = searchParams.getAll("sort") || []

  const filterParams = useMemo(() => {
    return {
      name: name,
      email: email,
    }
  }, [name, email])

  const sortParams = useMemo(() => {
    return {
      sort: sort,
    }
    //eslint-disable-next-line
  }, [JSON.stringify(sort)])

  const params = useMemo(() => {
    return {
      limit: Number(limit),
      offset: Number(offset),
    }
  }, [limit, offset])

  const [model, setModel] = useState<{
    filter: GridFilterModel
    sort: GridSortModel
  }>({
    filter: {
      items: [
        {
          field:
            Object.keys(filterParams).find(
              (key) => filterParams[key as keyof typeof filterParams],
            ) || "",
          operator: "contains",
          value: Object.values(filterParams).find((value) => value) || null,
        },
      ],
    },
    sort: [
      {
        field: sortParams.sort[0]?.replace("role", "role_id") || "",
        sort: sortParams.sort[1] as GridSortDirection,
      },
    ],
  })

  const { enqueueSnackbar } = useSnackbar()

  const handleClickVariant = (variant: VariantType, mess: string) => {
    enqueueSnackbar(mess, { variant })
  }

  useEffect(() => {
    if (!admin) navigate("/console")
    //eslint-disable-next-line
  }, [JSON.stringify(admin)])

  useEffect(() => {
    setModel({
      filter: {
        items: [
          {
            field:
              Object.keys(filterParams).find(
                (key) => filterParams[key as keyof typeof filterParams],
              ) || "",
            operator: "contains",
            value: Object.values(filterParams).find((value) => value) || null,
          },
        ],
      },
      sort: [
        {
          field: sortParams.sort[0]?.replace("role", "role_id") || "",
          sort: sortParams.sort[1] as GridSortDirection,
        },
      ],
    })
    //eslint-disable-next-line
  }, [sortParams, filterParams])

  useEffect(() => {
    if (newParams && newParams !== window.location.search.replace("?", "")) {
      setNewParams(window.location.search.replace("?", ""))
    }
    //eslint-disable-next-line
  }, [searchParams])

  useEffect(() => {
    if (newParams === window.location.search.replace("?", "")) return
    setParams(newParams)
    //eslint-disable-next-line
  }, [newParams])

  useEffect(() => {
    dispatch(getListUser({ ...filterParams, ...sortParams, ...params }))
    //eslint-disable-next-line
  }, [limit, offset, email, name, JSON.stringify(sort)])

  const getParamsData = () => {
    const dataFilter = Object.keys(filterParams)
      .filter((key) => filterParams[key as keyof typeof filterParams])
      .map((key) => `${key}=${filterParams[key as keyof typeof filterParams]}`)
      .join("&")
    return dataFilter
  }

  const paramsManager = useCallback(
    (page?: number) => {
      return `limit=${limit}&offset=${page ? page - 1 : offset}`
    },
    [limit, offset],
  )

  const handleSort = useCallback(
    (rowSelectionModel: GridSortModel) => {
      setModel({
        ...model,
        sort: rowSelectionModel,
      })
      let param
      const filter = getParamsData()
      if (!rowSelectionModel[0]) {
        param =
          filter || sortParams.sort[0] || offset
            ? `${filter ? `${filter}&` : ""}${paramsManager()}`
            : ""
      } else {
        param = `${filter}${
          rowSelectionModel[0]
            ? `${filter ? "&" : ""}sort=${rowSelectionModel[0].field.replace(
                "_id",
                "",
              )}&sort=${rowSelectionModel[0].sort}`
            : ""
        }&${paramsManager()}`
      }
      setNewParams(param)
    },
    //eslint-disable-next-line
    [paramsManager, getParamsData, model],
  )

  const handleFilter = (modelFilter: GridFilterModel) => {
    setModel({
      ...model,
      filter: modelFilter,
    })
    let filter = ""
    if (modelFilter.items[0]?.value) {
      filter = modelFilter.items
        .filter((item) => item.value)
        .map((item) => `${item.field}=${item?.value}`)
        .join("&")
    }
    const { sort } = sortParams
    const param =
      sort[0] || filter || offset
        ? `${filter}${
            sort[0] ? `${filter ? "&" : ""}sort=${sort[0]}&sort=${sort[1]}` : ""
          }&${paramsManager()}`
        : ""
    setNewParams(param)
  }

  const handleOpenModal = () => {
    setOpenModal(true)
  }

  type UserFormDTO = {
    id?: number
    name?: string
    email: string
    // temporarily use role's name (like "ADMIN") for select modal
    role_id?: string
  }
  const handleEdit = (dataEdit: UserFormDTO) => {
    setOpenModal(true)
    setDataEdit(dataEdit)
  }

  const onSubmitEdit = async (
    id: number | string | undefined,
    data: { [key: string]: string },
  ) => {
    const { role_id, ...newData } = data
    let newRole
    switch (role_id) {
      case "ADMIN":
        newRole = ROLE.ADMIN
        break
      case "OPERATOR":
        newRole = ROLE.OPERATOR
        break
    }
    if (id !== undefined) {
      const data = await dispatch(
        updateUser({
          id: id as number,
          data: { name: newData.name, email: newData.email, role_id: newRole },
          params: { ...filterParams, ...sortParams, ...params },
        }),
      )
      if (isRejectedWithValue(data)) {
        handleClickVariant("error", "Account update failed!")
        return
      } else {
        handleClickVariant(
          "success",
          "Your account has been edited successfully!",
        )
      }
    } else {
      const data = await dispatch(
        createUser({
          data: { ...newData, role_id: newRole } as AddUserDTO,
          params: { ...filterParams, ...sortParams, ...params },
        }),
      )
      if (isRejectedWithValue(data)) {
        if (!navigator.onLine) {
          handleClickVariant("error", "Account creation failed!")
          return
        }
        handleClickVariant("error", "This email already exists!")
        return
      } else {
        handleClickVariant(
          "success",
          "Your account has been created successfully!",
        )
      }
    }
    return
  }

  const handleOpenPopupDel = (id?: number, name?: string) => {
    if (!id) return
    setOpenDel({ id: id, name: name, open: true })
  }

  const handleClosePopupDel = () => {
    setOpenDel({ ...openDel, open: false })
  }

  const handleOkDel = async () => {
    if (!openDel?.id || !openDel) return
    const data = await dispatch(
      deleteUser({
        id: openDel.id,
        params: { ...filterParams, ...sortParams, ...params },
      }),
    )
    if (isRejectedWithValue(data)) {
      handleClickVariant("error", "Delete user failed!")
    } else {
      handleClickVariant("success", "Account deleted successfully!")
    }
    setOpenDel({ ...openDel, open: false })
  }

  const handleLimit = (event: ChangeEvent<HTMLSelectElement>) => {
    let filter = ""
    filter = Object.keys(filterParams)
      .filter((key) => filterParams[key as keyof typeof filterParams])
      .map((key) => `${key}=${filterParams[key as keyof typeof filterParams]}`)
      .join("&")
    const { sort } = sortParams
    const param = `${filter}${
      sort[0] ? `${filter ? "&" : ""}sort=${sort[0]}&sort=${sort[1]}` : ""
    }&limit=${Number(event.target.value)}&offset=0`
    setNewParams(param)
  }

  const handlePage = (event: ChangeEvent<unknown>, page: number) => {
    if (!listUser) return
    let filter = ""
    filter = Object.keys(filterParams)
      .filter((key) => filterParams[key as keyof typeof filterParams])
      .map((key) => `${key}=${filterParams[key as keyof typeof filterParams]}`)
      .join("&")
    const { sort } = sortParams
    const param = `${filter}${
      sort[0] ? `${filter ? "&" : ""}sort=${sort[0]}&sort=${sort[1]}` : ""
    }&limit=${listUser.limit}&offset=${(page - 1) * Number(limit)}`
    setNewParams(param)
  }

  const columns: GridColDef[] = [
    {
      headerName: "ID",
      field: "id",
      filterable: false,
      minWidth: 100,
      flex: 1,
    },
    {
      headerName: "Name",
      field: "name",
      type: "string",
      minWidth: 100,
      flex: 2,
      filterOperators: [
        {
          label: "Contains",
          value: "contains",
          getApplyFilterFn: () => null,
          InputComponent: ({ applyValue, item }: GridFilterInputValueProps) => {
            return (
              <Input
                sx={{ paddingTop: "16px" }}
                defaultValue={item.value || ""}
                onChange={(e) => {
                  if (timeout) clearTimeout(timeout)
                  timeout = setTimeout(() => {
                    applyValue({ ...item, value: e.target.value })
                  }, 300)
                }}
              />
            )
          },
        },
      ],
    },
    {
      headerName: "Role",
      field: "role_id",
      filterable: false,
      minWidth: 100,
      flex: 1,
      renderCell: (params: GridRenderCellParams) => {
        let role = ""
        switch (params.value) {
          case ROLE.ADMIN:
            role = "Admin"
            break
          case ROLE.OPERATOR:
            role = "OPERATOR"
            break
        }
        return <span>{role}</span>
      },
    },
    {
      headerName: "Mail",
      field: "email",
      type: "string",
      minWidth: 100,
      flex: 2,
      filterOperators: [
        {
          label: "Contains",
          value: "contains",
          getApplyFilterFn: () => null,
          InputComponent: ({ applyValue, item }: GridFilterInputValueProps) => {
            return (
              <Input
                sx={{ paddingTop: "16px" }}
                defaultValue={item.value || ""}
                onChange={(e) => {
                  if (timeout) clearTimeout(timeout)
                  timeout = setTimeout(() => {
                    applyValue({ ...item, value: e.target.value })
                  }, 300)
                }}
              />
            )
          },
        },
      ],
    },
    {
      headerName: "",
      field: "action",
      sortable: false,
      filterable: false,
      minWidth: 100,
      flex: 1,
      renderCell: (params: { row: UserDTO }) => {
        const { id, role_id, name, email } = params.row
        if (!id || !role_id || !name || !email) return null
        let role: string
        switch (role_id) {
          case ROLE.ADMIN:
            role = "ADMIN"
            break
          case ROLE.OPERATOR:
            role = "OPERATOR"
            break
        }

        return (
          <>
            <IconButton
              sx={{ color: "red" }}
              onClick={() =>
                handleEdit({ id, role_id: role, name, email } as UserFormDTO)
              }
            >
              <EditIcon sx={{ color: "black" }} />
            </IconButton>
            {!(params.row?.id === user?.id) ? (
              <IconButton
                sx={{ ml: 1.25 }}
                onClick={() =>
                  handleOpenPopupDel(params.row?.id, params.row?.name)
                }
              >
                <DeleteIcon sx={{ color: "red" }} />
              </IconButton>
            ) : null}
          </>
        )
      },
    },
  ]

  return (
    <AccountManagerWrapper>
      <Box
        sx={{
          display: "flex",
          justifyContent: "flex-end",
          gap: 2,
          marginBottom: 2,
        }}
      >
        <Button
          sx={{
            background: "#000000c4",
            "&:hover": { backgroundColor: "#00000090" },
          }}
          variant="contained"
          onClick={handleOpenModal}
        >
          Add
        </Button>
      </Box>
      <DataGrid
        sx={{ minHeight: 400, height: "calc(100vh - 300px)" }}
        columns={columns}
        rows={listUser?.items || []}
        filterMode={"server"}
        sortingMode={"server"}
        hideFooter
        onSortModelChange={handleSort}
        filterModel={model.filter}
        sortModel={model.sort as GridSortItem[]}
        onFilterModelChange={handleFilter}
      />
      {listUser && listUser.items.length > 0 ? (
        <PaginationCustom
          data={listUser}
          handlePage={handlePage}
          handleLimit={handleLimit}
          limit={Number(limit)}
        />
      ) : null}
      <PopupDelete
        open={openDel?.open || false}
        handleClose={handleClosePopupDel}
        handleOkDel={handleOkDel}
        name={openDel?.name}
      />
      {openModal ? (
        <ModalComponent
          onSubmitEdit={onSubmitEdit}
          setOpenModal={(flag) => {
            setOpenModal(flag)
            if (!flag) {
              setDataEdit({})
            }
          }}
          dataEdit={dataEdit}
        />
      ) : null}
      {loading ? <Loading /> : null}
    </AccountManagerWrapper>
  )
}

const AccountManagerWrapper = styled(Box)(({ theme }) => ({
  width: "80%",
  margin: theme.spacing(5, "auto"),
}))

const Modal = styled(Box)(() => ({
  position: "fixed",
  top: 0,
  left: 0,
  width: "100%",
  height: "100vh",
  display: "flex",
  justifyContent: "center",
  alignItems: "center",
  backgroundColor: "#cccccc80",
}))

const ModalBox = styled(Box)(() => ({
  width: 800,
  backgroundColor: "white",
  border: "1px solid black",
}))

const TitleModal = styled(Box)(({ theme }) => ({
  fontSize: 25,
  margin: theme.spacing(5),
}))

const BoxData = styled(Box)(() => ({
  marginTop: 35,
}))

const LabelModal = styled(Box)(({ theme }) => ({
  width: 300,
  display: "inline-block",
  textAlign: "end",
  marginRight: theme.spacing(0.5),
}))

const ButtonModal = styled(Box)(({ theme }) => ({
  button: {
    fontSize: 20,
  },
  display: "flex",
  justifyContent: "end",
  margin: theme.spacing(5),
}))

export default AccountManager
