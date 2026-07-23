一种基于网络增强的基因功能相似性计算方法（酵母数据）

项目描述
本代码库包含论文《一种基于网络增强的基因功能相似性计算方法》在酵母（Yeast）数据上的应用实现。方法基于酵母基因相互作用网络上的随机游走算法，结合基因本体（GO）注释与生化通路信息，计算基因功能相似性。

数据集信息
 原始数据文件
- `biochemical_pathways.tab` – 用于通路分析的酶学委员会（EC）分组数据（来源：http://sgd-archive.yeastgenome.org/）。
- `sgd.gaf` – 酵母基因的基因本体注释文件（来源：https://current.geneontology.org/products/pages/downloads.html）。
- `YeastNet.v3.txt` – 原始酵母基因相互作用网络（用于构建 `netHashMap`，来源：http://functionalnet.marcottelab.org/yeastnet/）。
- `YeastNet.v3.benchmark.txt` – 金标准基准基因对（用于生成 `geneList.txt`，来源：http://functionalnet.marcottelab.org/yeastnet/）。
- `go-basic.obo` – 基因本体术语总库文件（来源：https://www.geneontology.org/docs/download-ontology/）。

 中间数据文件（预处理过程中生成）
- `lfcPair.txt` – 用于 LFC 评分的基因对，来源于 `biochemical_pathways.tab` 的数据预处理。
- `geneList.txt` – 分析所需基因名称列表，从 `YeastNet.v3.benchmark.txt` 筛选得到。
- `netMatrix.txt` – 原始基因相互作用网络的矩阵表示。
- `pMatrix.txt` – 由 `RandWalk` 类计算得到的归一化概率转移矩阵。
- `child.txt` – GO 中后代术语的集合。
- `son.txt` – GO 中直接后代术语的集合。

 结果数据文件
- `result0.9.mat` – 重启参数为 0.9 时的随机游走结果矩阵。
- `similarityResult0.9.txt` – 重启参数为 0.9 时基因对的相似度得分。
- `lfc0.9.txt` – 各 EC 分组下的算法得分计算结果。

代码信息
本项目采用 Java 语言编写，主要类及其功能如下：
| 类名 | 功能描述 |
|---------|---------------|
| `Annotation` | GO 注释信息的数据模型。 |
| `FunctionNet` | 本体术语网络的数据模型。 |
| `LfcPair` | 获取 LFC 基因对的工具类。 |
| `Term` | GO 术语的数据模型。 |
| `Reader` | 特殊格式文件读取的辅助工具类。 |
| `Matrix` | 矩阵运算辅助工具类。 |
| `RandWalk` | 重启随机游走算法的实现类。 |
| `Calculator` | 核心计算类，按算法流程划分为四个部分，并提供三个主要外部接口。 |
| `Main` | 算法启动类，串联三个接口以运行整个流程。 |
| `ResultCalculate` | 结果数据测算类，用于初步分析结果特征。 |

使用说明
1. 导入项目
   打开 IntelliJ IDEA，选择`打开项目`。
   选择包含源代码的根目录，按默认选项完成导入。
   确保项目 SDK 设置为 JDK 8 或更高版本。

2. 配置数据路径
   确保上述“原始数据文件”中所列的全部文件放置于项目根目录下的 buf / 文件夹中
   （或按照代码中指定的路径存放）。其中，go-basic.obo直接置于该目录下，其他数据文件放置于`buf/Yeast/` 目录下。若路径不一致，请在 Main.java 或相关类中修改文件读取路径。

3. 选择术语领域
  注意Main.main()中的Calculator.data_initializer(gene_type, alpha, "molecular_function");
一行的"molecular_function"参数，这说明此时该程序将选择"molecular_function"作为术语领域， 可将其改为"biological_process"以更换术语领域。

4. 运行程序
   在项目视图中找到 Main 类（通常位于 src/ 目录下）。
   右键点击 Main 类，选择 Run 'Main.main()'。
   程序将依次执行以下步骤：
       加载 GO 结构与注释文件
       构建基因相互作用网络
       计算概率转移矩阵
       执行重启随机游走 (默认重启参数为 0.9)
       输出相似度结果文件与 LFC 评估文件

5. 结果输出
   结果文件将生成在results/文件夹或代码指定的输出目录下。

运行环境要求
IntelliJ IDEA
Java Development Kit (JDK) 8 或更高版本
矩阵运算需足够内存（完整酵母网络建议堆内存大于 4 GB）

方法步骤概述
算法流程在 Calculator 类中被划分为四个主要阶段：

数据加载与解析：读取 GO 层次结构（go-basic.obo）、注释文件（sgd.gaf）和相互作用网络（YeastNet.v3.txt）。

网络预处理：构建邻接矩阵并将其归一化为概率转移矩阵（pMatrix.txt）。

重启随机游走：利用 RandWalk 类对每一对基因执行随机游走计算。

结果评估：将相似度得分与基于通路的 EC 分组（LFC 评分）进行对比评估。

关于详细的数据预处理步骤，请参见论文正文“材料与方法”部分。
