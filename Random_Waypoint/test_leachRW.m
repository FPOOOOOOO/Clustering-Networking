function test_leachRW(s_mobility,s_input)
v_t = 0:s_input.TIME_STEP:s_input.SIMULATION_TIME;

%计算遮挡函数
    function isIntersect = checkLineIntersectSquare(x1,y1,x2,y2)
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
title('LEACH');
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
times = 1;                             %同一个mobility的实验次数
V_ele = 3*10^8;                        %光速

% 创建节点结构体、簇头结构体、记录数据结构体
S=repmat(struct('xd',0,'yd',0,'G',0','type','N','E',0,'ENERGY',0,'min_dis',0,'min_dis_cluster',0,'C',0,'T_Change',0),n+1,1);
C=repmat(struct('xd',0,'yd',0,'distance',0,'id',0,'child_num',0,'T_Change',0),n,1);% 记录当前轮簇头的相关数据
S(n+1).xd = sink.x;
S(n+1).yd = sink.y;
STATISTICS = repmat(struct('DEAD',0,'PACKETS_TO_BS',0,'PACKETS_TO_CH',0,'PACKETS_OUT_CH',0,'PACKETS_TO_NODES',0,'PACKETS_OUT_NODES',0,'CHNUMS',0,'T_Running',0,'T_Change',0,'E_Use',0),r,1);

% 这里是画图用的变量
FIRST_DEAD=zeros(1,times);% 记录每次实验第一次死掉节点时的轮数
%%%%%%% LEACH begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time = 1 : times
    %新实验次数需要重新归零的变量
    for i = 1:1:n                          %每个点的坐标及其初始化
        S(i).G = 0;                        % G<=0,说明该节点在候选集中
        S(i).type = 'N';                   % 普通节点
        S(i).E = Eo;                       % 初始能量
        S(i).ENERGY = 0;
        S(i).C = i;                        %初始簇头是自己
        C(i).xd = 0;                       %C(i)表示第i个簇头
        C(i).yd = 0;
        C(i).distance = 0;
        C(i).id = 0;                       %表示第i个簇头在S序列里的顺序
        C(i).child_num = 0;
    end
    t_HB_period = 0;
    for timeIndex = 1:1:length(v_t)
        t_HB_period = t_HB_period + 1
        %这个循环里就是轮次,显示时间
        t = v_t(timeIndex);
        set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
        %每个点的坐标及其初始化
        for i = 1:1:n
            if(S(i).E > 0)
                S(i).xd = vs_node(i).v_x(timeIndex);
                S(i).yd = vs_node(i).v_y(timeIndex);
            end
        end
        %遍历节点——查看死活
        for i=1:1:n
            if(S(i).E <= 0)                    %画出死亡节点
                set(vh_node_pos(i),'XData',S(i).xd,'YData',S(i).yd,'Marker','.','Color','red');
                STATISTICS(timeIndex).DEAD = STATISTICS(timeIndex).DEAD + 1;
            end

            if(S(i).E > 0)                     %画出没有死亡的节点
                S(i).type = 'N';
                set(vh_node_pos(i),'XData',S(i).xd,'YData',S(i).yd,'Marker','o','Color','blue');
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

        %LEACH的选举逻辑是t_last_period 200ms一次
        if(mod(t_HB_period,10) == 0) %10ms一次，进行.G的归零，和1/p一样，所以没问题
            %准备计算簇头数量
            for i = 1:1:n
                %重新选举数据清洗
                S(i).G = 0;
                S(i).type = 'N';
                S(i).C = i;
                C(i).xd = 0;
                C(i).yd = 0;
                C(i).distance = 0;
                C(i).id = 0;
                C(i).child_num = 0;
                %节点遍历,处理簇头数据
                if(S(i).E > 0)                     %节点还有能量
                    temp_rand = rand;              %随机种子
                    if((S(i).G) <= 0)
                        if(temp_rand <= (p/(1-p*mod(timeIndex,round(1/p)))))  %选择成为簇头
                            STATISTICS(timeIndex).CHNUMS = STATISTICS(timeIndex).CHNUMS + 1;
                            cluster = STATISTICS(timeIndex).CHNUMS;
                            S(i).type = 'C';
                            S(i).G = round(1/p) - 1;
                            C(cluster).xd = S(i).xd;
                            C(cluster).yd = S(i).yd;
                            set(vh_node_pos(i),'XData',S(i).xd,'YData',S(i).yd,'Marker','*','Color','green');
                            distance = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);  %簇头与BS的距离
                            C(cluster).distance = distance;
                            C(cluster).id = i;                                                   %第cluster个簇头点在S里的编号是i

                            %%簇头往各个节点广播包 CRC
                            distanceBroad = sqrt(xm*xm + ym*ym);
                            if(distanceBroad >= do)
                                S(i).E = S(i).E - (ETX * ctrPacketLength + Emp*ctrPacketLength*distanceBroad*distanceBroad*distanceBroad*distanceBroad);
                            else
                                S(i).E = S(i).E - (ETX * ctrPacketLength + Efc*ctrPacketLength*distanceBroad*distanceBroad);
                            end
                            STATISTICS(timeIndex).PACKETS_OUT_CH = STATISTICS(timeIndex).PACKETS_OUT_CH + 1;
                            %STATISTICS(timeIndex).PACKETS_TO_NODES = STATISTICS(timeIndex).PACKETS_TO_NODES + n - 1;                %这个NODES的能量还没算

                            %簇头向BS发送数据 DATA
                            distance;
                            if(distance >= do)
                                S(i).E = S(i).E - ((ETX+EDA)*packetLength + Emp*packetLength*(distance*distance*distance*distance ));
                            else
                                S(i).E = S(i).E - ((ETX+EDA)*packetLength + Efs*packetLength*(distance*distance));
                            end
                            STATISTICS(timeIndex).PACKETS_OUT_CH = STATISTICS(timeIndex).PACKETS_OUT_CH + 1;
                            STATISTICS(timeIndex).PACKETS_TO_BS = STATISTICS(timeIndex).PACKETS_TO_BS + 1;

                        end
                    end
                end
            end

            %节点遍历，处理非簇头数据
            for i = 1:1:n
                if(S(i).type == 'N' && S(i).E > 0)
                    min_dis = INFINITY;
                    if(STATISTICS(timeIndex).CHNUMS > 0)         %如果簇头存在
                        % 生成符合正态分布概率的随机数
                        mu = 75;  % 正态分布的均值
                        sigma = 10;  % 正态分布的标准差

                        min_dis_cluster = 1;    %距离最近的簇头编号
                        for c = 1:1:STATISTICS(timeIndex).CHNUMS   %遍历簇头找最近的
                            temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                            if(temp < min_dis)
                                min_dis = temp;
                                min_dis_cluster = c;
                            end
                            S(i).E = S(i).E - ERX * ctrPacketLength; %这里收下刚刚簇头广播的消息同时判断距离。 CRC
                            randomNumber = sigma * randn() + mu;
                            blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,C(c).xd,C(c).yd));
                            S(i).T_Change = S(i).T_Change + temp/V_ele * 1000 + blk * randomNumber;  %N：T_Change换算成ms
                            STATISTICS(timeIndex).PACKETS_TO_NODES = STATISTICS(timeIndex).PACKETS_TO_NODES + 1;
                        end

                        S(i).C = min_dis_cluster;  %标记自己连接的簇头是哪一个
                        randomNumber = sigma * randn() + mu;
                        blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,C(min_dis_cluster).xd,C(min_dis_cluster).yd));
                        C(min_dis_cluster).T_Change = C(min_dis_cluster).T_Change + 3 * min_dis/V_ele * 1000 + blk * randomNumber;  %簇头对其所连的的3K的 T_Change
                        C(min_dis_cluster).child_num = C(min_dis_cluster).child_num + 1;  %目前连接的簇头的子节点+1

                        %Energy dissipated by associated Cluster
                        %Head普通节点发送数据包到簇头消耗,和加入消息  CRC + DATA
                        if (min_dis > do)
                            S(i).E = S(i).E - (ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis)); %向簇头发送加入控制消息
                            S(i).E = S(i).E - (ETX*(packetLength) + Emp*packetLength*( min_dis * min_dis * min_dis * min_dis)); %向簇头数据包
                        else
                            S(i).E = S(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis * min_dis)); %向簇头发送加入控制消息
                            S(i).E = S(i).E -(ETX*(packetLength) + Efs*packetLength*( min_dis * min_dis)); %向簇头数据包
                        end
                        S(i).E = S(i).E - ERX*(ctrPacketLength);  %接收簇头确认加入控制消息 CRC
                        STATISTICS(timeIndex).PACKETS_OUT_NODES = STATISTICS(timeIndex).PACKETS_OUT_NODES + 2;
                        STATISTICS(timeIndex).PACKETS_TO_NODES = STATISTICS(timeIndex).PACKETS_TO_NODES + 1;
                        %Energy dissipated %簇头接收簇成员数据包消耗能量,接收加入消息和和确认加入消息
                        if(min_dis > 0) %这里是防止普通节点和簇头重合，但是会必然执行一遍 CRC + DATA
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA)*packetLength ); %接受簇成员发来的数据包
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ERX *ctrPacketLength ; %接收加入消息
                            if (min_dis > do)%簇头向簇成员发送确认加入的消息  CRC
                                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis));
                            else
                                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis * min_dis));
                            end
                            STATISTICS(timeIndex).PACKETS_TO_CH = STATISTICS(timeIndex).PACKETS_TO_CH + 2;
                            STATISTICS(timeIndex).PACKETS_OUT_CH = STATISTICS(timeIndex).PACKETS_OUT_CH + 1;
                        end
                        S(i).min_dis = min_dis;
                        S(i).min_dis_cluster = min_dis_cluster;
                    end
                end
            end

            %计算t_change & e_use
            t_change_N = 0;
            t_change_C = 0;
            e_use = 0;
            for i = 1:1:n
                t_change_N = max(t_change_N,S(i).T_Change);
                t_change_C = max(t_change_C,C(i).T_Change);
                e_use = e_use + (0.5 - S(i).E);
                S(i).T_Change = 0; %选举完就归零
                C(i).T_Change = 0;
            end
            STATISTICS(timeIndex).T_Change = t_change_N + t_change_C;
            STATISTICS(timeIndex).E_Use = e_use;
        end

        %准备处理非选举时刻的延时t_running
        for i = 1 : 1 : n
            if(S(i).type == 'N' && S(i).E > 0)
                % 生成符合正态分布概率的随机数
                mu = 75;  % 正态分布的均值
                sigma = 10;  % 正态分布的标准差
                randomNumber = sigma * randn() + mu;
                di = sqrt((C(S(i).C).xd - S(i).xd)^2 + (C(S(i).C).yd - S(i).yd)^2);
                blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,C(S(i).C).xd,C(S(i).C).yd));
                STATISTICS(timeIndex).T_Running = STATISTICS(timeIndex).T_Running + di/V_ele * 1000 +  blk * randomNumber; %以ms为单位
            end
            if(S(i).type == 'C' && S(i).E > 0)
                % 生成符合正态分布概率的随机数
                mu = 75;  % 正态分布的均值
                sigma = 10;  % 正态分布的标准差
                randomNumber = sigma * randn() + mu;
                di = sqrt((S(n+1).xd - S(i).xd)^2 + (S(n+1).yd - S(i).yd)^2);
                blk = (checkLineIntersectSquare(S(i).xd,S(i).yd,S(n+1).xd,S(n+1).yd));
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


    %将变量放到工作区
    assignin('base', 'STATISTICS', STATISTICS);
    assignin('base', 'S', S);
    assignin('base', 'C', C);

    %准备画图
    leach_x=1:1:timeIndex;
    leach_y=1:1:timeIndex;
    leach_z=1:1:timeIndex;
    leach_e=1:1:timeIndex;

    for i=1:1:timeIndex
        leach_x(i)= i;
        leach_y(i) = STATISTICS(i).T_Running;
        %z(i) = STATISTICS(i).PACKETS_TO_CH;
        leach_z(i) = STATISTICS(i).T_Change;
        leach_e(i) = STATISTICS(i).E_Use;
    end
    save('./LEACHR.mat', 'STATISTICS');
    assignin('base', 'leach_y', leach_y);
    assignin('base', 'leach_z', leach_z);
    assignin('base', 'leach_e', leach_e);
    %plot(x,y,'r',x,z,'b');
    figure(2);
    plot(leach_x,leach_y,'r');
    figure(3);
    plot(leach_x,leach_z,'b');
    figure(4);
    plot(leach_x,leach_e,'g');
    hold on;
end

end