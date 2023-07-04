import React, { FC, useState } from 'react'

import { GradientPickerPopover } from 'react-linear-gradient-picker'
import { SketchPicker, SketchPickerProps } from 'react-color'
import { PALETTE_COLOR_SHAPE_TYPE } from 'react-linear-gradient-picker'
import { ColorType } from 'store/slice/VisualizeItem/VisualizeItemType'

export const GradientColorPicker = React.memo<{
  colors: { rgb: string; offset: string }[]
  dispatchSetColor: (colorCode: ColorType[]) => void
}>(({ colors, dispatchSetColor }) => {
  const palette: PALETTE_COLOR_SHAPE_TYPE[] = colors.map((value) => {
    return {
      offset: value.offset,
      color: value.rgb,
    }
  })

  const onPaletteChange = (palette: PALETTE_COLOR_SHAPE_TYPE[]) => {
    const colorCode = palette.map((value) => {
      const color = value.color
      const rgbStr = color.replace(/[^0-9,]/g, '').split(',')
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
    <GradientPickerPopover
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
    </GradientPickerPopover>
  )
})

type WrapperPropTypes = {
  onSelect?: (color: string, opacity?: number) => void
  color?: string
}

const SketchPickerComponent: FC<SketchPickerProps> =
  SketchPicker as unknown as FC<SketchPickerProps>

const WrappedSketchPicker: React.FC<WrapperPropTypes> = ({
  onSelect,
  color,
}) => {
  return (
    <SketchPickerComponent
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
