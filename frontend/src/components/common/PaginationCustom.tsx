import {
  Box,
  FormControl,
  NativeSelect,
  Pagination,
  styled,
} from '@mui/material'
import { WorkspaceDataDTO } from '../../store/slice/Workspace/WorkspaceType'
import { ChangeEvent } from 'react'
import { UserListDTO } from '../../api/users/UsersApiDTO'

type PagiProps = {
  data: WorkspaceDataDTO | UserListDTO
  handlePage: (e: ChangeEvent<unknown>, page: number) => void
  handleLimit: (e: ChangeEvent<HTMLSelectElement>) => void
  limit: number | null
}

const PaginationCustom = ({
  data,
  handlePage,
  handleLimit,
  limit,
}: PagiProps) => {
  return (
    <PaginationCustomWrapper>
      <span>Rows per page: </span>
      <FormControl sx={{ width: 'max-content', margin: '4px 16px 0 12px' }}>
        <NativeSelect
          value={limit || 50}
          onChange={handleLimit}
          inputProps={{
            name: 'limit',
            id: 'uncontrolled-native',
          }}
        >
          <option value={10}>10</option>
          <option value={50}>50</option>
          <option value={100}>100</option>
        </NativeSelect>
      </FormControl>
      <span>{`${(data?.offset || 0) + 1} - ${
        (data?.offset || 0) + (data?.items?.length || 0)
      } of ${data?.total || 0}`}</span>
      <Pagination
        count={Math.ceil(data.total / data.limit)}
        page={Math.ceil(data.offset / data.limit) + 1}
        onChange={handlePage}
      />
    </PaginationCustomWrapper>
  )
}

const PaginationCustomWrapper = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'end',
  alignItems: 'center',
  marginTop: theme.spacing(2),
  gap: 2,
}))

export default PaginationCustom
