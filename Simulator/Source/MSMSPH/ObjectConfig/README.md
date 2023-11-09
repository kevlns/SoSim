1、数据组织格式：
``` json
{
  "unified_partRadius": 0.05,       // 求解器的统一粒子半径
  "objs": [                         // 物体列表，每个物体的基本属性相同
    {
      "mat": "fluid",           // 物体材质，具体选项见下
      "default": {              // 默认选项，在没有指定外部模型时，使用内部接口生成物体
        "shape": "cube",        // 物体形状
        "lb": [                 // left-bottom 坐标
          -1.5,
          -1.5,
          -1.5
        ],
        "size": [               // 物体尺寸
          3,
          3,
          3
        ]
      },
      "source_file": "",        // 外部模型路径，预支持ply模型。若有该项，则物体从模型加载而不使用 default 配置
      "m_partRadius": 0.05,     // 该物体的自有粒子半径
      "density_m": 1000,          // 物体的 rest-density_m
      "phase": 1,               // 物体的相序号（多相流相关）
      "velStart": [             // 物体初始速度
        0,
        0,
        0
      ]
    },
    {
     ...another obj...
    }
}
```
单个obj的属性值可以为空，但不能没有。


2、某些属性的选项如下：
``` json
    mat:    // 物体材质
        "fluid"
        "rigid"
        "elastic"
        "bound"

    shape:  // 默认选项中的物体形状
        "box"
        "cube"
        "plane-x"
        "cylinder"   // not yet
```