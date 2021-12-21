import { IJsonModel } from 'flexlayout-react'

export const PARAM_FORM_TABSET_ID = 'PARAM_FORM_TABSET'
export const DISPLAY_DATA_TABSET_ID = 'DISPLAY_DATA_TABSET'

const flexjson: IJsonModel = {
  global: {
    tabSetEnableDeleteWhenEmpty: false,
  },
  layout: {
    type: 'row',
    weight: 100,
    children: [
      {
        type: 'row',
        children: [
          {
            type: 'tabset',
            height: 300,
            selected: 0,
            children: [
              {
                type: 'tab',
                name: 'flowchart',
                component: 'flowchart',
                enableClose: false,
              },
            ],
          },
          {
            type: 'row',
            children: [
              {
                type: 'tabset',
                id: PARAM_FORM_TABSET_ID,
                selected: 0,
                enableMaximize: false,
                children: [],
              },
              {
                type: 'tabset',
                id: DISPLAY_DATA_TABSET_ID,
                selected: 0,
                width: 600,
                children: [],
              },
            ],
          },
        ],
      },
    ],
  },
  borders: [
    {
      type: 'border',
      location: 'left',
      size: 200,
      selected: 0,
      children: [
        {
          type: 'tab',
          name: 'sidebar',
          component: 'sidebar',
          enableClose: false,
        },
      ],
    },
  ],
}

export { flexjson }
