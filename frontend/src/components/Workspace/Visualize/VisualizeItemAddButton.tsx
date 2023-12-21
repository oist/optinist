import { FC } from "react"
import { useDispatch } from "react-redux"

import AddIcon from "@mui/icons-material/Add"
import Box from "@mui/material/Box"
import Button from "@mui/material/Button"
import Paper from "@mui/material/Paper"
import { styled } from "@mui/material/styles"

import {
  insertInitialItemToNextColumn,
  pushInitialItemToNewRow,
} from "store/slice/VisualizeItem/VisualizeItemSlice"

export const VisualizeItemAddButton: FC<{
  itemId?: number
}> = ({ itemId }) => {
  const dispatch = useDispatch()
  const onClick = () => {
    itemId != null
      ? dispatch(insertInitialItemToNextColumn(itemId))
      : dispatch(pushInitialItemToNewRow())
  }
  return (
    <StyledPaper elevation={0} variant="outlined">
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        height="100%"
      >
        <StyledButton onClick={onClick}>
          <AddIcon fontSize="large" color="primary" />
        </StyledButton>
      </Box>
    </StyledPaper>
  )
}

const StyledPaper = styled(Paper)(({ theme }) => ({
  width: 260,
  height: 255,
  border: "dashed",
  borderWidth: 2,
  borderColor: theme.palette.divider,
  margin: theme.spacing(1),
}))

const StyledButton = styled(Button)({
  width: "100%",
  height: "100%",
})
