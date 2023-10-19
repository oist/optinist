import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Input,
  styled,
} from "@mui/material"
import {
  DataGrid,
  GridRenderCellParams,
  GridValidRowModel,
} from "@mui/x-data-grid"
import {
  ChangeEvent,
  MouseEvent as MouseEventReact,
  useCallback,
  useEffect,
  useRef,
  useState,
} from "react"
import { useDispatch, useSelector } from "react-redux"
import CancelIcon from "@mui/icons-material/Cancel"
import { postListUserShareWorkspaces } from "store/slice/Workspace/WorkspaceActions"
import { selectListSearch, selectLoading } from "store/slice/User/UserSelector"
import { getListSearch } from "store/slice/User/UserActions"
import Loading from "components/common/Loading"
import { UserDTO } from "api/users/UsersApiDTO"
import CheckIcon from "@mui/icons-material/Check"
import { resetUserSearch } from "store/slice/User/UserSlice"
import { AppDispatch } from "../../store/store"

type PopupType = {
  open: boolean
  id: number
  handleClose: (v: boolean) => void
  title?: string
  usersShare?: {
    share_type?: number
    users: UserDTO[]
  }
}

type TableSearch = {
  usersSuggest: UserDTO[]
  onClose: () => void
  handleAddListUser: (user: UserDTO) => void
  stateUserShare: UserDTO[]
}

const TableListSearch = ({
  usersSuggest,
  onClose,
  handleAddListUser,
  stateUserShare,
}: TableSearch) => {
  const ref = useRef<HTMLLIElement | null>(null)

  useEffect(() => {
    window.addEventListener("mousedown", onMouseDown)
    return () => {
      window.removeEventListener("mousedown", onMouseDown)
    }
    //eslint-disable-next-line
  }, [])

  const onMouseDown = (event: MouseEvent) => {
    if (
      ref.current?.contains((event as any).target) ||
      (event as any).target.id === "inputSearch"
    )
      return
    onClose?.()
  }

  return (
    <TableListSearchWrapper ref={ref} onBlur={() => console.log(123)}>
      <UlCustom>
        {usersSuggest.map((item) => {
          const isSelected = stateUserShare.some((i) => i.id === item.id)
          return (
            <LiCustom
              key={item.id}
              onClick={() => handleAddListUser(item)}
              style={{
                cursor: isSelected ? "not-allowed" : "pointer",
              }}
            >
              {`${item.name} (${item.email})`}
              {isSelected ? <CheckIcon style={{ fontSize: 14 }} /> : null}
            </LiCustom>
          )
        })}
      </UlCustom>
    </TableListSearchWrapper>
  )
}
const PopupShare = ({
  open,
  handleClose,
  usersShare,
  id,
  title,
}: PopupType) => {
  const usersSuggest = useSelector(selectListSearch)
  const loading = useSelector(selectLoading)
  const [textSearch, setTextSearch] = useState("")
  const [stateUserShare, setStateUserShare] = useState(usersShare || undefined)
  const dispatch = useDispatch<AppDispatch>()
  let timeout = useRef<NodeJS.Timeout | undefined>()

  useEffect(() => {
    if (usersShare) {
      // setUserIdsSelected(usersShare.users.map(user => user.id));
      setStateUserShare(usersShare)
    }
  }, [usersShare])

  useEffect(() => {
    if (timeout.current) clearTimeout(timeout.current)
    if (!textSearch) {
      dispatch(resetUserSearch())
      return
    }
    timeout.current = setTimeout(() => {
      dispatch(getListSearch({ keyword: textSearch }))
    }, 300)
    //eslint-disable-next-line
  }, [textSearch])

  const handleShareFalse = (
    e: any,
    params: GridRenderCellParams<GridValidRowModel>,
  ) => {
    e.preventDefault()
    e.stopPropagation()
    if (!stateUserShare) return
    const indexCheck = stateUserShare.users.findIndex(
      (user) => user.id === params.id,
    )
    const newStateUserShare = stateUserShare.users.filter(
      (user, index) => index !== indexCheck,
    )
    setStateUserShare({ ...setStateUserShare, users: newStateUserShare })
  }

  const columnsShare = useCallback(
    (
      handleShareFalse: (
        e: MouseEventReact<HTMLButtonElement>,
        parmas: GridRenderCellParams<GridValidRowModel>,
      ) => void,
    ) => [
      {
        field: "name",
        headerName: "Name",
        minWidth: 140,
        renderCell: (params: GridRenderCellParams<GridValidRowModel>) => (
          <span>{params.row.name}</span>
        ),
      },
      {
        field: "email",
        headerName: "Email",
        minWidth: 280,
        renderCell: (params: GridRenderCellParams<GridValidRowModel>) => (
          <span>{params.row.email}</span>
        ),
      },
      {
        field: "share",
        headerName: "",
        filterable: false,
        sortable: false,
        minWidth: 130,
        renderCell: (params: GridRenderCellParams<GridValidRowModel>) => {
          if (!params.row.share) return ""
          return (
            <Button onClick={(e) => handleShareFalse(e, params)}>
              <CancelIcon color={"error"} />
            </Button>
          )
        },
      },
    ],
    //eslint-disable-next-line
    [JSON.stringify(stateUserShare?.users)],
  )

  const handleOke = async () => {
    if (!stateUserShare) return
    let newUserIds = stateUserShare.users.map((user) => user.id)
    await dispatch(
      postListUserShareWorkspaces({
        id,
        data: { user_ids: newUserIds as number[] },
      }),
    )
    handleClose(true)
  }

  const handleSearch = (event: ChangeEvent<HTMLInputElement>) => {
    setTextSearch(event.target.value)
  }

  const handleCloseSearch = () => {
    setTextSearch("")
    dispatch(resetUserSearch())
  }

  const handleAddListUser = (user: any) => {
    if (!usersSuggest || !stateUserShare) return
    if (!stateUserShare.users.find((item) => item.id === user.id)) {
      setStateUserShare({
        ...stateUserShare,
        users: [...stateUserShare.users, user],
      })
    }
  }

  const handleClosePopup = (event: any) => {
    if (event.key === "Escape") {
      handleClose(false)
    }
  }

  if (!usersShare) return null

  return (
    <Box>
      <DialogCustom
        open={open}
        onClose={handleClose}
        sx={{ margin: 0 }}
        onKeyDown={handleClosePopup}
      >
        <DialogTitle>{title || "Share Database Record"}</DialogTitle>
        <DialogContent>
          <>
            <Box style={{ position: "relative" }}>
              <Input
                id="inputSearch"
                sx={{ width: "60%" }}
                placeholder={"Search and add users"}
                value={textSearch}
                onChange={handleSearch}
              />
              {textSearch && usersSuggest ? (
                <TableListSearch
                  onClose={handleCloseSearch}
                  usersSuggest={usersSuggest}
                  stateUserShare={stateUserShare?.users || []}
                  handleAddListUser={handleAddListUser}
                />
              ) : null}
            </Box>
            <p>Permitted users</p>
            {stateUserShare && (
              <DataGrid
                sx={{ minHeight: 400 }}
                // onRowClick={handleShareTrue}
                rows={stateUserShare?.users.map((user: any) => ({
                  ...user,
                  share: true,
                }))}
                columns={columnsShare(handleShareFalse)}
                hideFooterPagination
              />
            )}
          </>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => handleClose(false)}>Cancel</Button>
          <Button onClick={handleOke}>Ok</Button>
        </DialogActions>
      </DialogCustom>
      {loading ? <Loading /> : null}
    </Box>
  )
}

const DialogCustom = styled(Dialog)(({ theme }) => ({
  "& .MuiDialog-container": {
    "& .MuiPaper-root": {
      width: "70%",
      maxWidth: "890px",
    },
  },
}))

const TableListSearchWrapper = styled(Box)(({ theme }) => ({
  position: "absolute",
  background: "#fff",
  zIndex: 100,
  width: "60%",
  boxShadow:
    "0 6px 16px 0 rgba(0,0,0,.08), 0 3px 6px -4px rgba(0,0,0,.12), 0 9px 28px 8px rgba(0,0,0,.05)",
  borderBottomLeftRadius: 8,
  borderBottomRightRadius: 8,
  maxHeight: 200,
  overflow: "auto",
}))

const UlCustom = styled("ul")(({ theme }) => ({
  listStyle: "none",
  padding: 0,
  margin: 0,
}))

const LiCustom = styled("li")(({ theme }) => ({
  padding: theme.spacing(1, 2),
  fontSize: 14,
  cursor: "pointer",
  display: "flex",
  justifyContent: "space-between",
  alignItems: "center",
  "&:hover": {
    backgroundColor: "rgba(0,0,0,.04)",
  },
}))

export default PopupShare
