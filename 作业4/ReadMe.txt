作业所用函数均存放于homework4.py中

ExampleFile中的.xlsx可移出与homework4.py放于同一路径，否则需要运行homework4_1_a()和homework4_1_b()后才可运行homework4_1_c()进行处理

其中homework4_1_a()函数生成作业第一题a问的样本与三角形，统计每个样本内出现不同三角形的数目，存放于H4a.xlsx中，总的三角形出现概率及其误差存放于Possibility-a.xlsx中；homework4_1_b()同理生成第一题b问中的样本与三角形，将统计结果存放于H4b.xlsx，三角形概率及误差存放于Possibility-b.xlsx内。

homework4_1_c()读取同一路径下的Possibility-a.xlsx和Possibility-b.xlsx，对其中的数据作比，调出比例偏离题目要求的三角形和其概率比例，输出到Ratio.xlsx文件中。

TriK()为用于挑选合适三角形数据输出为tex表格样式的脚本函数，可略。

homework4_2_a()输出蒙卡路径积分在[0,2]内均匀分布九个点的值，以及其与题中所给参考值对比图。

homework4_2_bc()生成第二题b问要求的组态，并处理成c问所要求的关联函数以及能量差。

Extra()计算10000个组态得到的能量差图样。