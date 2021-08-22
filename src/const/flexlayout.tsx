import { IJsonModel } from 'flexlayout-react'

var flexjson: IJsonModel = {
  global: {},
  layout: {
    type: 'row',
    weight: 100,
    children: [
      {
        type: 'row',
        weight: 10,
        children: [
          {
            type: 'tabset',
            selected: 0,
            children: [
              {
                type: 'tab',
                name: 'sidebar',
                component: 'sidebar',
              },
            ],
          },
          {
            type: 'tabset',
            selected: 0,
            children: [
              {
                type: 'tab',
                name: 'parameter',
                component: 'paramForm',
              },
            ],
          },
        ],
      },
      {
        type: 'tabset',
        weight: 50,
        selected: 0,
        children: [
          {
            type: 'tab',
            name: 'flow',
            component: 'flow',
          },
        ],
      },
    ],
  },
  borders: [
    {
      type: 'border',
      location: 'bottom',
      size: 100,
      children: [
        {
          type: 'tab',
          name: 'four',
          component: 'grid',
        },
      ],
    },
    {
      type: 'border',
      location: 'left',
      size: 100,
      children: [],
    },
  ],
}

export { flexjson }
