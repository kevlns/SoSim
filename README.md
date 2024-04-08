# SoSim

Physical Simulator of SoEngine

## Code Standards

1. 类名：采用 **UpperCamelCase** 格式，例如:

 ``` c++
    class HelloWorld {}
 ``` 

2. 变量名：采用**每单词分隔格式**，类中的非静态变量以 m_ 开头，静态变量直接命名，主机变量加 host_ 前缀，设备变量加 device_
   前缀，若有特别类型的设备变量，则加上特定的标识中缀，例如

```c++
    // member var
float m_name;
// static var
string path;
// host member var
Vec3f m_host_pos;
// cuda device var
float* m_device_cuda_age;
// opengl device var
Vec3f* m_device_gl_vertex;
```

3. 配置信息的结构体类型中的变量直接按照**每单词分隔格式**命名即可，例如

```c++
struct ObjectConfig {
    string name;
    string shape;
    float  life_time;
};
```

4. 文件命名时，按照**每单词分隔格式**并根据功能直观命名即可，例如 object_manager.hpp
5. 接口写注释，标明参数和返回结果的含义（ps：虽然俺也还没完全做到，但是新代码要注意哈）
6. 提交代码时，按照 **操作类型：简单描述** 的格式写提交信息，操作类型按照以下统一的格式

 ```c++
   // 新模块或者新功能
   feat: ...
   // 已有模块或功能上添加代码
   add: ...
   // 已有模块或功能上修改代码
   modify: ...
   // 在已推出的功能上修改bug
   fix: ...
   // 删除模块或功能的代码
   remove: ...
   // 合并分支
   merge: ...
```

## Build with CLion

### Windows

#### 1. setup

```bash
    git clone https://github.com/kevlns/SoSim.git
    git submodule init
    git submodule update
    cd /path/to/thirdparty/vcpkg
    .\bootstrap-vcpkg.bat
```

#### 2. Install dependencies

select vcpkg tool in the project  
![select_vcpkg.png](pics/select_vcpkg.png)  
vcpkg install  
![vcpkg_install.png](pics/vcpkg_install.png)

#### 3. Build and Run
