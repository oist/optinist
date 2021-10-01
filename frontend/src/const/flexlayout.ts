import { IJsonModel } from 'flexlayout-react'

const flexjson: IJsonModel = {
  global: {},
  layout: {
    type: 'row',
    weight: 100,
    children: [
      {
        type: 'row',
        children: [
          {
            type: 'tabset',
            height: 500,
            selected: 0,
            children: [
              {
                type: 'tab',
                name: 'flowchart',
                component: 'flowchart',
              },
            ],
          },
          {
            type: 'row',
            children: [
              {
                type: 'tabset',
                selected: 0,
                enableMaximize: false,
                children: [
                  {
                    type: 'tab',
                    name: 'parameter',
                    component: 'paramForm',
                  },
                ],
              },
              {
                type: 'tabset',
                selected: 0,
                enableDeleteWhenEmpty: false,
                children: [
                  {
                    type: 'tab',
                    name: 'image',
                    component: 'image',
                  },
                  {
                    type: 'tab',
                    name: 'output',
                    component: 'output',
                  },
                ],
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
