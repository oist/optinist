import { NWBListDTO, NWBListType } from './NWBType'

export function convertToNWBListType(dto: NWBListDTO) {
  const nwbList: NWBListType = {}
  Object.entries(dto).forEach(([name, value]) => {
    if (Object.prototype.hasOwnProperty.call(value, 'children')) {
      nwbList[name] = {
        type: 'parent',
        children: convertToNWBListType(
          (
            value as {
              children: NWBListDTO
            }
          ).children as NWBListDTO,
        ),
      }
    } else {
      nwbList[name] = {
        type: 'child',
        value,
      }
    }
  })
  return nwbList
}
