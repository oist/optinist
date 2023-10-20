import { FC, memo, useState } from "react"
import { SketchPicker } from "react-color"
import {
  GradientPickerPopover,
  PALETTE_COLOR_SHAPE_TYPE,
} from "react-linear-gradient-picker"

import { ColorType } from "store/slice/VisualizeItem/VisualizeItemType"

const GradientPickerPopoverComponent =
  // NOTE: use any because GradientPickerPopover type definition is not provided
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  GradientPickerPopover as unknown as FC<any>

interface GradientColorPickerProps {
  colors: { rgb: string; offset: string }[]
  dispatchSetColor: (colorCode: ColorType[]) => void
}

export const GradientColorPicker = memo(function GradientColorPicker({
  colors,
  dispatchSetColor,
}: GradientColorPickerProps) {
  const palette: PALETTE_COLOR_SHAPE_TYPE[] = colors.map((value) => {
    return {
      offset: value.offset,
      color: value.rgb,
    }
  })

  const onPaletteChange = (palette: PALETTE_COLOR_SHAPE_TYPE[]) => {
    const colorCode = palette.map((value) => {
      const color = value.color
      const rgbStr = color.replace(/[^0-9,]/g, "").split(",")
      const rgb = {
        r: Number(rgbStr[0]),
        g: Number(rgbStr[1]),
        b: Number(rgbStr[2]),
      }
      return {
        rgb: `rgb(${rgb.r}, ${rgb.g}, ${rgb.b})`,
        offset: value.offset,
      }
    })
    dispatchSetColor(colorCode)
  }

  const [open, setOpen] = useState(false)

  return (
    <GradientPickerPopoverComponent
      open={open}
      setOpen={() => setOpen(!open)}
      // showAnglePicker={true}
      width={150}
      maxStops={10}
      paletteHeight={25}
      palette={palette}
      onPaletteChange={onPaletteChange}
      flatStyle={true}
    >
      <WrappedSketchPicker />
    </GradientPickerPopoverComponent>
  )
})

type WrapperPropTypes = {
  onSelect?: (color: string, opacity?: number) => void
  color?: string
}

const WrappedSketchPicker: FC<WrapperPropTypes> = ({ onSelect, color }) => {
  return (
    <SketchPicker
      color={color}
      width="150px"
      // styles={{width: "10px"}}
      onChange={(c) => {
        const { r, g, b, a } = c.rgb
        onSelect?.(`rgb(${r}, ${g}, ${b})`, a)
      }}
    />
  )
}
