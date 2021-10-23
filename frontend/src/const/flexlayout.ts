import { IJsonModel } from 'flexlayout-react'
import { toLayoutTab } from 'utils/FlexLayoutUtils'
import {
  INITIAL_IMAGE_ELEMENT_ID,
  INITIAL_ALGO_ELEMENT_ID,
  INITIAL_ALGO_ELEMENT_NAME,
  INITIAL_IMAGE_ELEMENT_NAME,
} from './flowchart'

export const PARAM_FORM_TABSET_ID = 'PARAM_FORM_TABSET_ID'
export const OUTPUT_TABSET_ID = 'OUTPUT_TABSET_ID'

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
                id: PARAM_FORM_TABSET_ID,
                selected: 0,
                enableMaximize: false,
                children: [
                  toLayoutTab(
                    INITIAL_ALGO_ELEMENT_ID,
                    'paramForm',
                    INITIAL_ALGO_ELEMENT_NAME,
                  ),
                ],
              },
              {
                type: 'tabset',
                id: OUTPUT_TABSET_ID,
                selected: 0,
                enableDeleteWhenEmpty: false,
                children: [
                  toLayoutTab(
                    INITIAL_IMAGE_ELEMENT_ID,
                    'image',
                    INITIAL_IMAGE_ELEMENT_NAME,
                  ),
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
