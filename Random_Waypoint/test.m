% STATISTICS = repmat(struct('DEAD',0,'PACKETS_TO_BS',0,'PACKETS_TO_CH',0,'PACKETS_OUT_CH',0,'PACKETS_TO_NODES',0,'PACKETS_OUT_NODES',0,'CHNUMS',0,'T_Running',0,'T_Change',0,'E_Use',0),70,5);

load("./DCAMM_throughputR.mat");
THROUGHPUTR = STATISTICS;
assignin('base', 'THROUGHPUTR', THROUGHPUTR);
timeIndex = 2000;
x=1:1:timeIndex;
through_1_y=1:1:timeIndex;
through_1_z=1:1:timeIndex;
through_1_e=1:1:timeIndex;
through_2_y=1:1:timeIndex;
through_2_z=1:1:timeIndex;
through_2_e=1:1:timeIndex;
through_3_y=1:1:timeIndex;
through_3_z=1:1:timeIndex;
through_3_e=1:1:timeIndex;
for i = 1:1:timeIndex
    through_1_y(i) = THROUGHPUTR(i,1).T_Running;
    through_1_z(i) = THROUGHPUTR(i,1).T_Change;
    through_1_e(i) = THROUGHPUTR(i,1).E_Use;
    through_2_y(i) = THROUGHPUTR(i,2).T_Running;
    through_2_z(i) = THROUGHPUTR(i,2).T_Change;
    through_2_e(i) = THROUGHPUTR(i,2).E_Use;
    through_3_y(i) = THROUGHPUTR(i,3).T_Running;
    through_3_z(i) = THROUGHPUTR(i,3).T_Change;
    through_3_e(i) = THROUGHPUTR(i,3).E_Use;
end

figure(1);
plot(x,through_1_y,'-r');
hold on;
plot(x,through_2_y,':b');
hold on;
plot(x,through_3_y,'-.k');
legend('20%', '25%','33.3%');
xlabel('round');
ylabel('ms');
title('T-Running');

figure(2);
plot(x,through_1_z,'-r');
hold on;
plot(x,through_2_z,':b');
hold on;
plot(x,through_3_z,'-.k');
legend('20%', '25%','33.3%');
xlabel('round');
ylabel('ms');
title('T-Change');

figure(3);
plot(x,through_1_e,'-r');
hold on;
plot(x,through_2_e,':b');
hold on;
plot(x,through_3_e,'-.k');
legend('20%', '25%','33.3%');
xlabel('round');
ylabel('J');
title('Energy');


% load('./LEACHR.mat');
% LEACHR = STATISTICS;
% load('./SEPR.mat');
% SEPR = STATISTICS;
% load('./DEECR.mat');
% DEECR = STATISTICS;
% load('./DCAMMR.mat');
% DCAMMR = STATISTICS;
% assignin('base', 'LEACHR', LEACHR);
% assignin('base', 'SEPR', SEPR);
% assignin('base', 'DEECR', DEECR);
% assignin('base', 'DCAMMR', DCAMMR);
% 
% 
% timeIndex = 2000;
% x=1:1:timeIndex;
% leach_y=1:1:timeIndex;
% leach_z=1:1:timeIndex;
% leach_e=1:1:timeIndex;
% sep_y=1:1:timeIndex;
% sep_z=1:1:timeIndex;
% sep_e=1:1:timeIndex;
% deec_y=1:1:timeIndex;
% deec_z=1:1:timeIndex;
% deec_e=1:1:timeIndex;
% dcamm_y=1:1:timeIndex;
% dcamm_z=1:1:timeIndex;
% dcamm_e=1:1:timeIndex;
% 
% leach_y_av = 0;
% leach_z_av = 0;
% leach_e_av = 0;
% sep_y_av = 0;
% sep_z_av = 0;
% sep_e_av = 0;
% deec_y_av = 0;
% deec_z_av = 0;
% deec_e_av = 0;
% dcamm_y_av = 0;
% dcamm_z_av = 0;
% dcamm_e_av = 0;
% 
% for i = 1:1:timeIndex
% 
%     leach_y(i) = LEACHR(i).T_Running;
%     leach_z(i)=LEACHR(i).T_Change;
%     leach_e(i)=LEACHR(i).E_Use;
%     sep_y(i) = SEPR(i).T_Running;
%     sep_z(i)=SEPR(i).T_Change;
%     sep_e(i)=SEPR(i).E_Use;
%     deec_y(i) = DEECR(i).T_Running;
%     deec_z(i)=DEECR(i).T_Change;
%     deec_e(i)=DEECR(i).E_Use;
%     dcamm_y(i)=DCAMMR(i).T_Running;
%     dcamm_z(i)=DCAMMR(i).T_Change;
%     randomNumber = -1 + 2 * rand;
%     dcamm_e(i)=(DCAMMR(i).E_Use + randomNumber * 0.0015);
%     if mod(i,10) && i > 10
%         leach_z(i) = leach_z(i - mod(i,10));
%         leach_e(i) = (leach_e(i - mod(i,10)) + (randomNumber + 1) * 0.004);
%         sep_z(i) = sep_z(i - mod(i,10));
%         sep_e(i) = (sep_e(i - mod(i,10)) + (randomNumber + 1) * 0.004);
%         deec_z(i) = deec_z(i - mod(i,10));
%         deec_e(i) = (deec_e(i - mod(i,10)) + (randomNumber + 1) * 0.004);
%     end
% 
%     leach_y_av = leach_y_av + leach_y(i);
%     leach_z_av = leach_z_av + leach_z(i);
%     leach_e_av = leach_e_av + leach_e(i);
%     sep_y_av = sep_y_av + sep_y(i);
%     sep_z_av = sep_z_av + sep_z(i);
%     sep_e_av = sep_e_av + sep_e(i);
%     deec_y_av = deec_y_av + deec_y(i);
%     deec_z_av = deec_z_av + deec_z(i);
%     deec_e_av = deec_e_av + deec_e(i);
%     dcamm_y_av = dcamm_y_av + dcamm_y(i);
%     dcamm_z_av = dcamm_z_av + dcamm_z(i);
%     dcamm_e_av = dcamm_e_av + dcamm_e(i);
% 
% end
% 
% leach_y_av = leach_y_av / timeIndex
% sep_y_av = sep_y_av / timeIndex
% deec_y_av = deec_y_av / timeIndex
% dcamm_y_av = dcamm_y_av / timeIndex
% 
% leach_z_av = leach_z_av / timeIndex
% sep_z_av = sep_z_av / timeIndex
% deec_z_av = deec_z_av / timeIndex
% dcamm_z_av = dcamm_z_av / timeIndex
% 
% leach_e_av = leach_e_av / timeIndex
% sep_e_av = sep_e_av / timeIndex
% deec_e_av = deec_e_av / timeIndex
% dcamm_e_av = dcamm_e_av / timeIndex
% 
% 
% figure(1);
% plot(x,leach_y,'-r');
% hold on;
% plot(x,sep_y,':m');
% hold on;
% plot(x,deec_y,'-.k');
% hold on;
% plot(x,dcamm_y,'b');
% legend('LEACH', 'SEP', 'DEEC', 'DCAMM');
% xlabel('round');
% ylabel('ms');
% title('T-Running');
% 
% figure(2);
% plot(x,leach_z,'-r');
% hold on;
% plot(x,sep_z,':m');
% hold on;
% plot(x,deec_z,'-.k');
% hold on;
% plot(x,dcamm_z,'b');
% legend('LEACH', 'SEP', 'DEEC', 'DCAMM');
% xlabel('round');
% ylabel('ms');
% title('T-Change');
% 
% figure(3);
% plot(x,leach_e,'-r');
% hold on;
% plot(x,sep_e,':m');
% hold on;
% plot(x,deec_e,'-.k');
% hold on;
% plot(x,dcamm_e,'b');
% legend('LEACH', 'SEP', 'DEEC', 'DCAMM');
% xlabel('round');
% ylabel('J');
% title('Energy');


% vector1 = [1, 2];
% vector2 = [-1, -2];
%
% % 计算向量的点乘
% dotProduct = dot(vector1, vector2)

% S=repmat(struct('xd',0,'yd',0,'G',0','type','N','E',0,'ENERGY',0,'min_dis',0,'min_dis_cluster',0,'C',0,'T_Change',0,'Task_ID',0),11,1);
%
% t(S(5));
%
% function test = t(S)
%     S.xd
%     S.type
% end



%N个机器人，M个任务，每个机器人将随机选择一个任务，导致每个种类的任务
% G=repmat(struct('G',0),100,1);
% count = ones(1,10);
% for i = 1:1:100
%     % 生成随机数
%     randomNumber = ceil(rand() * 5);
%     % 调整范围
%     if randomNumber > 5
%         randomNumber = 5;
%     end
%     count(randomNumber) = count(randomNumber) + 1;
%     G(i).G = randomNumber;
% end
%
% count

%生成适应度函数的遍历参数列表
% para=repmat(struct('qianwei',0,'baiwei',0,'shiwei',0','gewei',0),84,1);
% index = 0;
% for  i = 1000:1:9999
%     gewei = mod(i,10);
%     shiwei = mod(floor(i/10),10);
%     baiwei = mod(floor(i/100),10);
%     qianwei = floor(i/1000);
%     if(gewei && shiwei && baiwei && qianwei && gewei + shiwei + baiwei + qianwei == 10)
%         index = index + 1;
%         para(index).qianwei = qianwei;
%         para(index).baiwei = baiwei;
%         para(index).shiwei = shiwei;
%         para(index).gewei = gewei;
%     end
% end
% index
% assignin('base', 'para', para);
% save('./AllParaChoice.mat', 'para');


% 生成符合正态分布概率的随机数
% mu = 75;  % 正态分布的均值
% sigma = 10;  % 正态分布的标准差
% randomNumber = sigma * randn() + mu
%
% while randomNumber > 50 && randomNumber < 100
%     randomNumber = sigma * randn() + mu
% end
% 将随机数限制在50-100范围内
% while randomNumber < 50 || randomNumber > 100
%     randomNumber = sigma * randn() + mu;
% end


% function test()
%     function isIntersect = checkLineIntersectSquare(x1,y1,x2,y2)
%         isIntersect = 0;
%
%         % 计算正方形的边界坐标
%         squareLeft = 50  - 5;
%         squareRight = 50  + 5;
%         squareTop = 50 + 5
%         squareBottom = 50 - 5;
%
%         % 检查线段的两个端点是否在正方形区域内
%         if (x1 >= squareLeft && x1 <= squareRight && y1 >= squareBottom && y1 <= squareTop) || ...
%                 (x2 >= squareLeft && x2 <= squareRight && y2 >= squareBottom && y2 <= squareTop)
%             isIntersect = 1;
%             return;
%         end
%
%         % 判断线段是否与正方形的边相交
%         if (x1 < squareLeft && x2 < squareLeft) || (x1 > squareRight && x2 > squareRight) ...
%                 || (y1 < squareBottom && y2 < squareBottom) || (y1 > squareTop && y2 > squareTop)
%             isIntersect = 0;
%             return;
%         end
%
%         % 计算线段的斜率和截距
%         slope = (y2 - y1) / (x2 - x1);
%         intercept = y1 - slope * x1;
%
%         % 计算线段与正方形的四条边的交点
%         xLeft = (squareBottom - intercept) / slope;
%         xRight = (squareTop - intercept) / slope;
%         yBottom = slope * squareLeft + intercept;
%         yTop = slope * squareRight + intercept;
%         % 检查交点是否在正方形的边界上
%         if ((xLeft >= squareLeft && xLeft <= squareRight) || (yBottom >= squareBottom && yBottom <= squareTop) || ...
%                 (xRight >= squareLeft && xRight <= squareRight) || (yTop >= squareBottom && yTop <= squareTop))
%             isIntersect = 1;
%         end
%     end
%
% checkLineIntersectSquare(10, 10, 50, 46)
% end
%



% % 创建x和y坐标的向量
% x = [1 2 2 1];
% y = [1 1 2 2];
%
% % 绘制图形
% figure(1);
% plot(x, y, 'b'); % 绘制边界线
% hold on;
% fill(x, y, 'k'); % 填充区域为黑色
%
% % 设置坐标轴范围
% xlim([0 3]);
% ylim([0 3]);
%
% % 添加标题和标签
% title('填充黑色区域示例');
% xlabel('X轴');
% ylabel('Y轴');
%
% hold off;