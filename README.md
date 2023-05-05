## 230407
1、动态分簇组网机制!搞这个！nnd！\
2、全视距+智能反射面\ 
3、速度高/切换速度快（时刻在切） \
4、通信质量函数来预判是否切换 \
5、参考手机的基站切换（磁悬浮 MCS 控制 BTS  \
6、raft、paxos \

## 230409
1、LEACH，与SEP、DEEC、SHC比较，这三个可以去看看\
2、分簇算法:\
LEACH、TEEN、PEGASIS、SEP、EEUC、FLEACH、SSEP、SHC、DeepClus、DEEC、HEAD、EEUC、UCFIA
LEACH（Low-Energy Adaptive Clustering Hierarchy）- 2000 \
TEEN（Threshold-sensitive Energy Efficient sensor Network protocol）-2002 \
PEGASIS（Power-Efficient Gathering in Sensor Information Systems）-2002 \
SEP（Stable Election Protocol）-2005 \
EEUC（Energy-Efficient Unequal Clustering） - 2006 \
FLEACH（Fuzzy LEACH）- 2006 \
SSEP（Self-Stabilizing Energy-efficient Protocol）-2007 \
DEEC（Distributed Energy-Efficient Clustering） - 2007 \
SHC（Spatially Homogeneous Clustering） - 2012 \
DeepClus（Deep Clustering for Unsupervised Learning in WSNs）-2019 \
HEED（Hybrid Energy-Efficient Distributed Clustering\
UCFIA（Unequal Clustering Algorithm for WSN based on Fuzzy Logic and Improved ACO）\

3、选中继、
4、群智能算法：
遗传算法（Genetic Algorithm，GA）- 1960 \
免疫算法（Immune Algorithm，IA）- 1970 \
蚁群算法（Ant Colony Optimization，ACO） - 1991 \
粒子群算法（Particle Swarm Optimization，PSO）- 1995 \
人工蜂群算法（Artificial Bee Colony，ABC） -2005 \

## 230410
1、ABC选簇头，然后LEACH加入簇头，簇内TDMA，簇间CDMA
2、Kmeans或者FCM分簇，然后选簇头
3、基于能量均衡移动多层分簇算法，Load_Balance Mutil_Driven clustering algorithm with mobile nodes.
4、Mutil_Layer ——Mobility-Aware Hierarchical Clustering in Mobile Wireless Sensor Networks + 能量平衡的无线传感器网络数据采集动态分簇算法
5、Energy_Balance —— 能量平衡的无线传感器网络数据采集动态分簇算法

## 230411
1、LEACH的matlab研究
2、FCM和Fuzzing Grid待看

## 230412
1、LEACH的matlab能跑了\
2、簇头选择——基于能量均衡高效的LEACH路由协议优化策略_黄利晓，或者根据我的K值来算

## 230413
1、算法总设计目标：低延时、减少冗余、能量总消耗少
2、基本假设：能量供应充足，延时最低（如何仿真延时）

## 230414
1、延时，20ms变完结构
2、用功能分簇，然后选簇头
3、遮挡就让他断开

## 230417
1、重点放到延迟身上，要把这个过程厘清
2、如果按照WiFi 802.11 DCF来的话就是CSMA/DA
3、如果按照之前的文章来的话，就是TDMA

## 230419
1、RWP跑起来了，要结合一下，它是先把所有点的轮次模拟好，然后在渲染，应该根据时间，边走边渲染

## 230426
1、继续阅读matlab代码\
2、LEACH好像有节点冲突，但是RWP里好像没有这个情况 —— 错了错了是有这个情况的，就是在移动所以没看太出来\
3、每次把随机出来的部分保存一份csv，这样有对比性\

## 230504
1、进行一个数据的保存\
2、R = unifrnd(A, B, m, n) 生成m * n 的数组\

## 230505
1、进行一个数据的保存\
2、
