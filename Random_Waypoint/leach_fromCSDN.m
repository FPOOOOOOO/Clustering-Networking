%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEP: A Stable Election Protocol for clustered                        %
%      heterogeneous wireless sensor networks                          %
% (c) Georgios Smaragdakis                                             %
% WING group, Computer Science Department, Boston University           %
% You can find full documentation and related information at:          %
% http://csr.bu.edu/sep                                                %
% To report your comment or any bug please send e-mail to:             %
% gsmaragd@cs.bu.edu                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the LEACH [1] code we have used.                             %
% The same code can be used for FAIR if m=1                            %
% [1] W.R.Heinzelman, A.P.Chandrakasan and H.Balakrishnan,             %
%     "An application-specific protocol architecture for wireless      %
%      microsensor networks"                                           %
%     IEEE Transactions on Wireless Communications, 1(4):660-670,2002  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm = 100;
ym = 100;
%x and y Coordinates of the Sink
%sink.x =0.5 * xm;
%sink.y = ym + 50;
sink.x=50;
sink.y=50;
%sink.x=0.5*xm;
%sink.y=0.5*ym;

%Number of Nodes in the field
n = 100
%Optimal Election Probability of a node to become cluster head
p=0.05;
packetLength =6400;%数据包长度
ctrPacketLength = 200;%控制包长度
%Energy Model (all values in Joules)
%Initial Energy
Eo = 0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

INFINITY = 999999999999999;
%maximum number of rounds
rmax=2001
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do

do=sqrt(Efs/Emp); %do全程没变过，是多径衰落的模型里的传输距离分界线，和能量有关

%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;%坐标 要从mobility里拿来用
    %XR(i)=S(i).xd;       % XR（i）是什么
    S(i).yd=rand(1,1)*ym;
    %YR(i)=S(i).yd;
    S(i).G=0;            % .G<=0, 说明该节点在候选集中
    %initially there are no cluster heads only nodes
    S(i).type='N';%普通节点
    S(i).E=Eo;
    S(i).ENERGY=0;
    plot(S(i).xd,S(i).yd,'o','color',[0.3 0.3 1]);
    hold on;
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x'); %画BS

%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax %主循环,每次1轮,这个后续要合并成每个时刻，因为这样才能保证运动的点进行簇头选择
    r;
    %Operation for epoch
    if(mod(r, round(1/p))==0)%类似epoch的概念，在p等于0.05的情况下是每20个r然后对所有的.G清零
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;  %.cl是什么
        end
    end

    hold off;

    %Number of dead nodes
    dead=0;

    %counter for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS=0;
    packets_TO_CH=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads per round
    %Upper Alphabet means per round（r）
    PACKETS_TO_CH(r+1)=0;
    PACKETS_TO_BS(r+1)=0;

    figure(1);

    for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            plot(S(i).xd,S(i).yd,'red .');
            dead=dead+1;
            hold on;
        end

        if (S(i).E>0)
            S(i).type='N';
            plot(S(i).xd,S(i).yd,'o','color',[0.3 0.3 1]);
            hold on;
        end
        %在这里加一个循环，应该可以得到每一轮完之后，网络剩余的能量
    end

    if (dead == n)%节点全部死亡退出循环
        break;
    end

    STATISTICS(r+1).DEAD=dead;
    DEAD(r+1)=dead;

    %When the first node dies
    if (dead==1)
        if(flag_first_dead==0)
            first_dead=r; %第一个节点死亡时的轮数
            flag_first_dead=1;
        end
    end

    countCHs=0;
    cluster=1; %第一个簇头
    for i=1:1:n  %对点遍历，当前轮次
        if(S(i).E>0)
            temp_rand=rand; %还有能量的节点里随机
            if ((S(i).G)<=0) %如果该节点在候选集合中
                %Election of Cluster Heads
                if(temp_rand <=(p/(1-p*mod(r,round(1/p)))))
                    countCHs = countCHs+1;

                    S(i).type = 'C';
                    S(i).G = round(1/p)-1;    %这不是一个固定值吗？为什么要减1
                    C(cluster).xd = S(i).xd;  %簇头的位置
                    C(cluster).yd = S(i).yd;
                    plot(S(i).xd,S(i).yd,'k*'); %黑色的节点

                    distance=sqrt((S(i).xd-(S(n+1).xd))^2 + (S(i).yd-(S(n+1).yd))^2);%到sink的距离

                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    %X(cluster)=S(i).xd;
                    %Y(cluster)=S(i).yd;
                    cluster=cluster+1;

                    %广播自成为簇头
                    distanceBroad = sqrt(xm*xm+ym*ym);   %对角线
                    if (distanceBroad >=do)
                        S(i).E = S(i).E-(ETX*ctrPacketLength + Emp*ctrPacketLength*(distanceBroad*distanceBroad*distanceBroad*distanceBroad));%广播自成为簇头
                    else
                        S(i).E = S(i).E-(ETX*ctrPacketLength + Efs*ctrPacketLength*(distanceBroad*distanceBroad));
                    end
                    %Calculation of Energy dissipated 簇头自己发送数据包能量消耗
                    distance;
                    if(distance>=do)
                        S(i).E = S(i).E-((ETX+EDA)*packetLength+ Emp*packetLength*(distance*distance*distance*distance ));
                    else
                        S(i).E = S(i).E-((ETX+EDA)*packetLength+ Efs*packetLength*(distance*distance));
                    end
                    packets_TO_BS = packets_TO_BS+1;
                    PACKETS_TO_BS(r+1) = packets_TO_BS;
                end
            end
        end
    end

    STATISTICS(r+1).CLUSTERHEADS = cluster-1;%统计第r轮簇头数目,r是从0开始的,所以加1;cluster最后要-1,是上面的循环多加了1
    CLUSTERHS(r+1)= cluster-1;

    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
        if (S(i).type=='N' && S(i).E>0) %普通节点
            % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%默认距离是到sink的距离
            min_dis = INFINITY;
            if(cluster-1>=1)%如果有簇头存在
                min_dis_cluster = 1;
                %加入最近的簇头
                for c = 1:1:cluster-1 %簇头数量一共是cluster-1，遍历簇头找最近的
                    %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
                    temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                    if (temp<min_dis)
                        min_dis = temp;
                        min_dis_cluster = c;
                    end
                    %接收簇头发来的广播的消耗
                    S(i).E = S(i).E - ETX * ctrPacketLength;
                end


                %Energy dissipated by associated Cluster Head普通节点发送数据包到簇头消耗,和加入消息
                min_dis;
                if (min_dis > do)
                    S(i).E = S(i).E - (ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis)); %向簇头发送加入控制消息
                    S(i).E = S(i).E - (ETX*(packetLength) + Emp*packetLength*( min_dis * min_dis * min_dis * min_dis)); %向簇头数据包
                else
                    S(i).E = S(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis * min_dis)); %向簇头发送加入控制消息
                    S(i).E = S(i).E -(ETX*(packetLength) + Efs*packetLength*( min_dis * min_dis)); %向簇头数据包
                end
                S(i).E = S(i).E - ETX*(ctrPacketLength);  %接收簇头确认加入控制消息

                %Energy dissipated %簇头接收簇成员数据包消耗能量,接收加入消息和和确认加入消息
                if(min_dis > 0)
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA)*packetLength ); %接受簇成员发来的数据包
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ERX *ctrPacketLength ; %接收加入消息
                    if (min_dis > do)%簇头向簇成员发送确认加入的消息
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis));
                    else
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis * min_dis));
                    end
                    PACKETS_TO_CH(r+1) = n - dead - cluster + 1; %所有的非死亡的普通节点都发送数据包
                end

                S(i).min_dis = min_dis;
                S(i).min_dis_cluster = min_dis_cluster;

            end
        end
    end
    %hold on;

    countCHs;
    rcountCHs = rcountCHs + countCHs;
end

x=1:1:r;
y=1:1:r;
%z=1:1:r;

figure(2);
for i=1:1:r;
    x(i)=i;
    y(i) = n - STATISTICS(i).DEAD;
    %z(i)=CLUSTERHS(i);
end
%plot(x,y,'r',x,z,'b');
plot(x,y,'r');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round
%  %通过这一个量，可以判断10%节点死亡时经过的时间
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died  %这是第一个节点死亡时经过的轮数
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%标准数据：sink(50,175) ,ctrPacketLength=200,packetLength=4000,Eo=2J.