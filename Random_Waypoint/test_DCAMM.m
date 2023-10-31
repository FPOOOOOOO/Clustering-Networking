function test_DCAMM(s_mobility,s_input)
v_t = 0:s_input.TIME_STEP:s_input.SIMULATION_TIME;

for nodeIndex = 1:s_mobility.NB_NODES
    %Simple interpolation (linear) to get the position, anytime.
    %Remember that "interp1" is the matlab function to use in order to
    %get nodes' position at any continuous time.
    %插值得到v_x,这个是每个时间点的vx,根据时间来的。s_mobility里的是根据逻辑来的。
    vs_node(nodeIndex).v_x = interp1(s_mobility.VS_NODE(nodeIndex).V_TIME,s_mobility.VS_NODE(nodeIndex).V_POSITION_X,v_t);
    vs_node(nodeIndex).v_y = interp1(s_mobility.VS_NODE(nodeIndex).V_TIME,s_mobility.VS_NODE(nodeIndex).V_POSITION_Y,v_t);
end

figure;
hold on;
%vh_nodes_pos 在工作区是Line，这里对每个节点进行操作，打印初始化位置，vh_nodes_pos是用来在time轴上进行变换的
for nodeIndex = 1:s_mobility.NB_NODES
    vh_node_pos(nodeIndex) = plot(vs_node(nodeIndex).v_x(1),vs_node(nodeIndex).v_y(1),'o','color',[0.3 0.3 1]);
end

%画遮挡的区域
for obstacle = 1:s_input.OBSTAClE_NUM
    half_len = s_input.OBSTACLE_EDGE / 2;
    tempx = [s_mobility.OBSTACLE(obstacle).V_POSITION_X-half_len s_mobility.OBSTACLE(obstacle).V_POSITION_X+half_len s_mobility.OBSTACLE(obstacle).V_POSITION_X+half_len s_mobility.OBSTACLE(obstacle).V_POSITION_X-half_len];
    tempy = [s_mobility.OBSTACLE(obstacle).V_POSITION_Y-half_len s_mobility.OBSTACLE(obstacle).V_POSITION_Y-half_len s_mobility.OBSTACLE(obstacle).V_POSITION_Y+half_len s_mobility.OBSTACLE(obstacle).V_POSITION_Y+half_len];
    plot(tempx,tempy,'b');
    hold on;
    fill(tempx,tempy,'k');
end

plot(50,50,'red x');         %画BS

title(cat(2,'Simulation time (sec): ',num2str(s_mobility.SIMULATION_TIME)));
xlabel('X (meters)');
ylabel('Y (meters)');
title('DCAMM');
%ht在工作区是Text,标注
ht = text(s_input.V_POSITION_X_INTERVAL(1),s_input.V_POSITION_Y_INTERVAL(2),cat(2,'Time (sec) = 0'));
%ht = text(min(vs_node(1).v_x),max(vs_node(1).v_y),cat(2,'Time (sec) = 0'));
%伸缩坐标轴
axis([s_input.V_POSITION_X_INTERVAL(1) s_input.V_POSITION_X_INTERVAL(2) s_input.V_POSITION_Y_INTERVAL(1) s_input.V_POSITION_Y_INTERVAL(2)]);
%axis([min(vs_node(1).v_x) max(vs_node(1).v_x) min(vs_node(1).v_y) max(vs_node(1).v_y)]);
hold off;

%%%%% LEACH paras %%%%%
xm = s_input.V_POSITION_X_INTERVAL(2); %场地大小X轴
ym = s_input.V_POSITION_Y_INTERVAL(2); %场地大小y轴
sink.x = 50;                           %网关节点x
sink.y = 50;                           %网关节点y
n = s_input.NB_NODES;                  %移动节点数量
p = 0.1;                               %LEACH的选取概率
packetLength = 6400;                   %数据包长度
ctrPacketLength = 200;                 %控制包长度
Eo = 0.5;                              %初始能量
ETX = 50*0.000000001;                  %Eelec=Etx=Erx
ERX = 50*0.000000001;                  %
Efs = 10*0.000000000001;               %Transmit Amplifier types
Emp = 0.0013*0.000000000001;           %
EDA = 5*0.000000001;                   %Data Aggregation Energy
INFINITY = 999999999999999;
r = length(v_t);                       %轮次从mobility里面拿
do=sqrt(Efs/Emp);                      %do全程没变过，是多径衰落的模型里的传输距离分界线，和能量有关
V_ele = 3*10^8;                        %光速
para_index = 84;                        %四个不同参数的遍历顺序，先从1开始
data = load('./AllParaChoice.mat');
para = data.para;
times = 1;                             %同一个mobility的实验次数

% 创建节点结构体、簇头结构体、记录数据结构体
S=repmat(struct('xd',0,'yd',0,'G',0','type','N','E',0,'min_dis',0,'min_dis_cluster',0,'C',0,'T_Change',0,'Task_ID',0,'P',-999),n+1,1);
C=repmat(struct('xd',0,'yd',0,'distance',0,'id',0,'child_num',0,'T_Change',0,'P',-999),n,1);% 记录当前轮簇头的相关数据
S(n+1).xd = sink.x;
S(n+1).yd = sink.y;
STATISTICS = repmat(struct('DEAD',0,'PACKETS_TO_BS',0,'PACKETS_TO_CH',0,'PACKETS_OUT_CH',0,'PACKETS_TO_NODES',0,'PACKETS_OUT_NODES',0,'CHNUMS',0,'T_Running',0,'T_Change',0,'E_Use',0),r,1);
t_change_Av = zeros(1,84);

% 这里是画图用的变量
FIRST_DEAD=zeros(1,times);% 记录每次实验第一次死掉节点时的轮数
%%%%%%% LEACH begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time = 1 : times
    %新实验次数需要重新归零的变量
    task_id_nums = zeros(1,s_input.TASK_NUM); %每种任务的节点个数
    for i = 1:1:n                          %每个点的坐标及其初始化
        S(i).G = 0;                        % G<=0,说明该节点在候选集中
        S(i).type = 'N';                   % 普通节点
        S(i).E = Eo;                       % 初始能量
        S(i).C = i;                        %初始簇头是自己
        S(i).T_change = 0;
        S(i).Task_ID = s_mobility.VS_NODE(i).TASK_ID; %当前节点的任务簇
        task_id_nums(S(i).Task_ID) = task_id_nums(S(i).Task_ID) + 1;  %每个任务簇内含有的节点数量
        S(i).P = 0;                        %一开始的簇头适应度函数值为0
        C(i).xd = 0;                       %C（i）表示Task_ID为i的任务簇的簇头
        C(i).yd = 0;
        C(i).distance = 0;
        C(i).id = 0;                       %Task_ID为i的任务簇的簇头在S中的序列号
        C(i).child_num = 0;
        C(i).T_Change = 0;
    end

    t_HB_period = 0;
    for timeIndex = 1:1:length(v_t)
        time
        t_HB_period = t_HB_period + 1
        %这个循环里就是轮次,显示时间
        t = v_t(timeIndex);
        %set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
        %每个点的坐标及其初始化
        for i = 1:1:n
            if(S(i).E > 0)
                S(i).xd = vs_node(i).v_x(timeIndex);
                S(i).yd = vs_node(i).v_y(timeIndex);
            end
        end
        %遍历节点——查看死活,同时计算活着节点的适应度函数
        for i=1:1:n
            if(S(i).E <= 0)                    %画出死亡节点
                %set(vh_node_pos(i),'XData',S(i).xd,'YData',S(i).yd,'Marker','.','Color','red');
                STATISTICS(timeIndex).DEAD = STATISTICS(timeIndex).DEAD + 1;
            end

            if(S(i).E > 0)                     %画出没有死亡的节点
                S(i).type = 'N';
                S(i).P = PiCalculate(S,i,para,para_index,s_mobility,s_input);
                %set(vh_node_pos(i),'XData',S(i).xd,'YData',S(i).yd,'Marker','o','Color','blue');
            end
        end
        %全死了就结束
        if(STATISTICS(timeIndex).DEAD == n)
            break;
        end
        %如果有节点死了，那么就记录第一个死亡的节点时的轮次
        if(STATISTICS(timeIndex).DEAD && ~FIRST_DEAD(time))
            FIRST_DEAD(time) = timeIndex;
        end

        %节点遍历，广播簇头适应度函数并收下其他簇内节点的广播包
        for i = 1:1:n
            now_task = S(i).Task_ID;
            %节点遍历,处理簇头数据
            if(S(i).E > 0)                     %节点还有能量
                if S(i).P > C(now_task).P
                    C(now_task).P = S(i).P;
                    C(now_task).id = i;
                end
                %%节点算出来以后往各个节点广播包 CRC，是发送包%
                distanceBroad = sqrt(xm*xm + ym*ym);
                if(distanceBroad >= do)
                    S(i).E = S(i).E - (ETX * ctrPacketLength + Emp*ctrPacketLength*distanceBroad*distanceBroad*distanceBroad*distanceBroad);
                else
                    S(i).E = S(i).E - (ETX * ctrPacketLength + Efc*ctrPacketLength*distanceBroad*distanceBroad);
                end
                STATISTICS(timeIndex).PACKETS_OUT_NODES = STATISTICS(timeIndex).PACKETS_OUT_NODES + 1;                        %这时候还不知道是簇头还是节点，所以没法标记，全部算成节点，取消C的包统计
                %%当前节点要收下其他Gi - 1 + 1个的包;
                S(i).E = S(i).E - ERX * ctrPacketLength * task_id_nums(now_task);
            end
        end

        %%标记簇头的位置，并在S序列中标记其的type
        for i = 1:1:s_input.TASK_NUM
            if C(i).id > 0   %说明有这个任务簇的节点
                C(i).xd = S(C(i).id).xd;
                C(i).yd = S(C(i).id).yd;
                S(C(i).id).type = 'C';
            end
        end

        %节点遍历，和簇头连接，簇头广播，并计算t_change,目前最后的那个1没有算
        for i = 1 : 1 : n
            if(S(i).E > 0)
                % 生成符合正态分布概率的随机数
                mu = 75;  % 正态分布的均值
                sigma = 10;  % 正态分布的标准差
                randomNumber = sigma * randn() + mu;
                di = sqrt((C(S(i).Task_ID).xd - S(i).xd)^2 + (C(S(i).Task_ID).yd - S(i).yd)^2);
                blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,C(S(i).Task_ID).xd,C(S(i).Task_ID).yd,s_mobility,s_input));
                S(i).T_Change = S(i).T_Change + di/V_ele * 1000 + blk * randomNumber;

            end
        end

        %计算t_change & E_use
        t_change_N = 0;
        t_change_C = 0;
        e_use = 0;
        for i = 1:1:n
            t_change_N = max(t_change_N,S(i).T_Change);
            S(i).T_Change = 0; %选举完就归零
            C(i).T_Change = 0;
            e_use = e_use + (0.5 - S(i).E);
        end
        t_change_Av(time) = t_change_Av(time) + t_change_N;
        STATISTICS(timeIndex).T_Change = t_change_N + t_change_C;
        STATISTICS(timeIndex).E_Use = e_use;


        %准备处理非选举时刻的延时t_running
        for i = 1 : 1 : n
            if(S(i).type == 'N' && S(i).E > 0)
                % 生成符合正态分布概率的随机数
                mu = 75;  % 正态分布的均值
                sigma = 10;  % 正态分布的标准差
                randomNumber = sigma * randn() + mu;
                di = sqrt((C(S(i).Task_ID).xd - S(i).xd)^2 + (C(S(i).Task_ID).yd - S(i).yd)^2);
                blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,C(S(i).Task_ID).xd,C(S(i).Task_ID).yd,s_mobility,s_input));
                STATISTICS(timeIndex).T_Running = STATISTICS(timeIndex).T_Running + di/V_ele * 1000 +  blk * randomNumber; %以ms为单位
            end
            if(S(i).type == 'C' && S(i).E > 0)
                % 生成符合正态分布概率的随机数
                mu = 75;  % 正态分布的均值
                sigma = 10;  % 正态分布的标准差
                randomNumber = sigma * randn() + mu;
                di = sqrt((S(n+1).xd - S(i).xd)^2 + (S(n+1).yd - S(i).yd)^2);
                blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,S(n+1).xd,S(n+1).yd,s_mobility,s_input));
                STATISTICS(timeIndex).T_Running = STATISTICS(timeIndex).T_Running + di/V_ele * 1000 +  blk * randomNumber; %以ms为单位
            end
        end

        %%%%%%% LEACH end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         for nodeIndex = 1:s_mobility.NB_NODESg
        %             set(vh_node_pos(nodeIndex),'XData',vs_node(nodeIndex).v_x(timeIndex),'YData',vs_node(nodeIndex).v_y(timeIndex));
        %         end
        %当前轮次结束，更新figure(1)内容
        drawnow;
    end

    t_change_Av(time) = t_change_Av(time) / timeIndex;

    %将变量放到工作区
    assignin('base', 'STATISTICS', STATISTICS);
    assignin('base', 'S', S);
    assignin('base', 'C', C);

    %准备画图
    dcamm_x=1:1:timeIndex;
    dcamm_y=1:1:timeIndex;
    dcamm_z=1:1:timeIndex;
    dcamm_e=1:1:timeIndex;


    for i=1:1:timeIndex
        dcamm_x(i)= i;
        dcamm_y(i) = STATISTICS(i).T_Running;  %节点存活数
        dcamm_z(i) = STATISTICS(i).T_Change;
        dcamm_e(i) = STATISTICS(i).E_Use;
    end
    save('./DCAMMR.mat', 'STATISTICS');
    figure(2);
    plot(dcamm_x,dcamm_y,'r');
    figure(3);
    plot(dcamm_x,dcamm_z,'b');
    figure(4);
    plot(dcamm_x,dcamm_e,'g');
    hold on;
end
%
% %画图找最好的84分支1
% x=1:1:time;
% y=1:1:time;
% ymin = 9999;
% yminindex = 0;
% figure(3);
% for i = 1:1:time
%     x(i) = i;
%     y(i) = t_change_Av(i);
%     if y(i) < ymin
%         ymin = y(i);
%         yminindex = i;
%     end
% end
% yminindex
% plot(x,y,'g');
% hold on;
end

%计算遮挡函数
function isIntersect = checkLineIntersectSquare(x1,y1,x2,y2,s_mobility,s_input)
isIntersect = 0;
for obi = 1: s_input.OBSTAClE_NUM
    % 计算正方形的边界坐标
    squareLeft = s_mobility.OBSTACLE(obi).V_POSITION_X  - s_mobility.OBSTACLE_EDGE /2;
    squareRight = s_mobility.OBSTACLE(obi).V_POSITION_X  + s_mobility.OBSTACLE_EDGE /2;
    squareTop = s_mobility.OBSTACLE(obi).V_POSITION_Y  + s_mobility.OBSTACLE_EDGE /2;
    squareBottom = s_mobility.OBSTACLE(obi).V_POSITION_Y  - s_mobility.OBSTACLE_EDGE /2;


    % 检查线段的两个端点是否在正方形区域内
    if (x1 >= squareLeft && x1 <= squareRight && y1 >= squareBottom && y1 <= squareTop) || ...
            (x2 >= squareLeft && x2 <= squareRight && y2 >= squareBottom && y2 <= squareTop)
        isIntersect = 1;
        return;
    end

    % 判断线段是否与正方形的边相交
    if (x1 < squareLeft && x2 < squareLeft) || (x1 > squareRight && x2 > squareRight) ...
            || (y1 < squareBottom && y2 < squareBottom) || (y1 > squareTop && y2 > squareTop)
        isIntersect = 0;
        continue;
    end

    % 计算线段的斜率和截距
    slope = (y2 - y1) / (x2 - x1);
    intercept = y1 - slope * x1;

    % 计算线段与正方形的四条边的交点
    xLeft = (squareBottom - intercept) / slope;
    xRight = (squareTop - intercept) / slope;
    yBottom = slope * squareLeft + intercept;
    yTop = slope * squareRight + intercept;
    % 检查交点是否在正方形的边界上
    if ((xLeft >= squareLeft && xLeft <= squareRight) || (yBottom >= squareBottom && yBottom <= squareTop) || ...
            (xRight >= squareLeft && xRight <= squareRight) || (yTop >= squareBottom && yTop <= squareTop))
        isIntersect = 1;
        return;
    end
end
end

%计算适应度函数,这里输入的应该是S和其对应的i，para和其对应的j
function PiValue = PiCalculate(S,i,para,j,s_mobility,s_input)
n = s_input.NB_NODES;
now_task = S(i).Task_ID;
%计算di mi rssi mri
alpha1 = para(j).qianwei;
alpha2 = para(j).baiwei;
alpha3 = para(j).shiwei;
alpha4 = para(j).gewei;

pdi = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);

pcenter_xd = 0;
pcenter_yd = 0;
now_task_count = 0;
for k = 1:1:n
    if(S(k).Task_ID == now_task)
        now_task_count = now_task_count + 1;
        pcenter_xd = pcenter_xd + S(k).xd;
        pcenter_yd = pcenter_yd + S(k).yd;
    end
end
pcenter_xd = pcenter_xd / now_task_count;
pcenter_yd = pcenter_yd / now_task_count;
%这里的pmi的速度向量比较难算，就用和簇节点加权重心的差值来算
pmi = sqrt((S(i).xd - pcenter_xd)^2 + (S(i).yd - pcenter_yd)^2);

%     其计算公式为
% d=10^((ABS(RSSI)-A)/(10*n))
% 其中
% A为距离设备1米时的rssi 绝对值。
% n为环境衰减因子。
% Python的代码实现为

% def rssiRange(rssi):
% A = -38.0
% n = 27.0
% iRssi = abs(rssi)
% power = float(( iRssi - A ) / ( 10 * n ))
% result = pow(10 , power)
% return result
pA = -31;
pshuaijian = 27;
pblk = checkLineIntersectSquare(S(i).xd,S(i).yd,S(n+1).xd,S(n+1).yd,s_mobility,s_input);
prssi = -(log10(pdi) * 10 * pshuaijian  + pA) - pblk * 55;
vector_now = [S(i).xd, S(i).yd];
vector_next = [S(i+1).xd, S(i+1).yd];
L_now = pdi;
L_next = sqrt((S(i + 1).xd - S(n+1).xd)^2 + (S(i + 1).yd - S(n+1).yd)^2);

pmri = acos(dot(vector_now,vector_next) / L_now / L_next);

PiValue = alpha1 / pdi + alpha2 /pmi + alpha3 * prssi + alpha4 / pmri;
end