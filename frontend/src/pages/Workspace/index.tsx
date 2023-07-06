// import { useEffect } from 'react'
import { useSelector /*, useDispatch */ } from 'react-redux'
import { Box, styled } from '@mui/material'
import {
  DataGrid,
  GridColDef,
  GridRenderCellParams,
} from '@mui/x-data-grid'
import { Link } from 'react-router-dom'
import Loading from '../../components/common/Loading'
import {
  selectIsLoadingWorkspaceList,
  selectWorkspaceList,
} from 'store/slice/Workspace/WorkspaceSelector'

const columns: GridColDef[] = [
  {
    field: 'workspace_id',
    headerName: 'ID',
    renderCell: (params: GridRenderCellParams<string>) => (
      <Link to={`/workspaces/${params.value}`}>{params.value}</Link>
    ),
  },
]

const Workspaces = () => {
  // const dispatch = useDispatch()
  const workspaces = useSelector(selectWorkspaceList)
  const loading = useSelector(selectIsLoadingWorkspaceList)

  /* TODO: Add get workspace apis and actions
  useEffect(() => {
    dispatch(getWorkspaceList())
    //eslint-disable-next-line
  }, [])
  */

  return (
    <WorkspacesWrapper>
      <WorkspacesTitle>Workspaces</WorkspacesTitle>
      <DataGrid
        autoHeight
        checkboxSelection
        rows={workspaces.map((ws) => ({
          id: ws.workspace_id,
          workspace_id: ws.workspace_id,
        }))}
        columns={columns}
      />
      {loading ? <Loading /> : null}
    </WorkspacesWrapper>
  )
}

const WorkspacesWrapper = styled(Box)(({ theme }) => ({
  width: '100%',
  padding: theme.spacing(2),
  height: 'calc(100% - 90px)',
  overflow: 'auto',
}))

const WorkspacesTitle = styled('h1')(({ theme }) => ({}))

export default Workspaces
