import { TextField, TextFieldProps } from "@mui/material"
import { styled } from "@mui/material/styles"

export const ParamTextField = styled((props: TextFieldProps) => (
  <TextField variant="outlined" fullWidth {...props} />
))(({ theme }) => ({
  "& .MuiInputLabel-root": {
    "&.MuiInputLabel-formControl": {
      position: "static",
      transform: "none",
      transition: "none",
      fontWeight: "bold",
      fontSize: "0.95rem",
      color: theme.palette.text.secondary,
    },
  },
  "& .MuiOutlinedInput-root": {
    marginTop: 0,
  },
  "& .MuiOutlinedInput-input": {
    paddingTop: "10px",
    paddingBottom: "8px",
    height: "auto",
  },
  "& .MuiOutlinedInput-notchedOutline": {
    top: 0,
    "& legend": {
      display: "none",
    },
  },
}))
