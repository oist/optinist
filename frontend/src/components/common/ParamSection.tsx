import { memo, ReactNode } from "react"

import { InsertDriveFileOutlined } from "@mui/icons-material"
import {
  Box,
  Chip,
  Divider,
  InputLabel,
  styled,
  Tooltip,
  Typography,
  TypographyProps,
} from "@mui/material"

import { getLabelByPath } from "store/slice/FlowElement/FlowElementUtils"

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

interface FileNameChipProps {
  filePath: string | null
}

export const FileNameChip = memo(function FileNameChip({
  filePath,
}: FileNameChipProps) {
  const label = filePath ? getLabelByPath(filePath) : "No file is selected"
  return (
    <Tooltip title={label}>
      <Chip
        icon={<InsertDriveFileOutlined />}
        label={label}
        sx={{ marginBottom: 2 }}
      />
    </Tooltip>
  )
})
