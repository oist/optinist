declare module 'react-linear-gradient-picker' {
  import React from 'react'
  export declare type PALETTE_COLOR_SHAPE_TYPE = {
    id?: number
    color: string
    offset: string
    opacity?: number
  }
  declare type GRADIENT_PICKER_PROP = {
    onPaletteChange: (value: PALETTE_COLOR_SHAPE_TYPE[]) => void
    paletteHeight?: number
    width?: number
    stopRemovalDrop?: number
    maxStops?: number
    minStops?: number
    flatStyle?: boolean
    palette?: PALETTE_COLOR_SHAPE_TYPE[]
  }
  declare type ANGLE_PICKER_PROP = {
    angle?: number
    setAngle?: (angle: number) => void
    size?: number
    snap?: number
  }
  declare type GRADIENT_PICKER_POPOVER_PROP = GRADIENT_PICKER_PROP &
    ANGLE_PICKER_PROP & {
      showAnglePicker?: boolean
      open: boolean
      setOpen: (open: boolean) => void
      trigger?: Function
    }
  export declare const GradientPicker: React.FC<GRADIENT_PICKER_PROP>
  export declare const GradientPickerPopover: React.FC<GRADIENT_PICKER_POPOVER_PROP>
}
