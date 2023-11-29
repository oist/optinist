import { memo, ReactNode } from "react"

import {
  Box,
  Divider,
  InputLabel,
  Typography,
  TypographyProps,
} from "@mui/material"
import { styled } from "@mui/material/styles"

interface SectionProps {
  title: string
  children: ReactNode
}

export const ParamSection = memo(function Section({
  title,
  children,
}: SectionProps) {
  return (
    <Box marginBottom={2}>
      <SectionTitle>{title}</SectionTitle>
      {children}
      <Divider style={{ marginTop: 8 }} />
    </Box>
  )
})

export const SectionTitle = styled((props: TypographyProps) => (
  <Typography {...props} />
))(({ theme }) => ({
  fontWeight: "bold",
  fontSize: "1.1rem",
  marginBottom: theme.spacing(1),
}))

export const FieldLabel = styled(InputLabel)({
  fontWeight: "bold",
  fontSize: "0.95rem",
  color: "text.secondary",
})
