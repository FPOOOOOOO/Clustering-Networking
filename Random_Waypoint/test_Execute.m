%Testing Random Waypoint mobility model.
clear all;clc;close all;

s_input = struct('V_POSITION_X_INTERVAL',[0 100],...    %(m)
                 'V_POSITION_Y_INTERVAL',[0 100],...    %(m)
                 'V_SPEED_INTERVAL',[0.2 22],...        %(m/s)
                 'V_PAUSE_INTERVAL',[0 7],...           %pause time (s)
                 'V_WALK_INTERVAL',[2.00 6.00],...      %walk time (s)
                 'V_DIRECTION_INTERVAL',[-180 180],...  %(degrees)
                 'SIMULATION_TIME',20,...               %(s)
                 'NB_NODES',20,...                      %节点数
                 'OBSTAClE_NUM',4,...                   %障碍物数量
                 'OBSTACLE_EDGE', 10,...                %长方形障碍物的边长
                 'TASK_NUM',4,...                       %任务数量
                 'TIME_STEP',0.01);                      %
% s_mobility = Generate_Mobility(s_input);
% % 生成带有系统时间的文件名
% currentDateTime = datetime('now', 'Format', 'MMddHHmm');
% currentDateTimeStr = datestr(currentDateTime, 'mmddHHMM');
% fileName = ['./mobility/mobility_', currentDateTimeStr, '.mat'];
% % 将结构体变量保存到带有系统时间的文件名中
% save(fileName, 's_mobility');

load("mobility/mobility_09262052.mat");

%test_Animate(s_mobility,s_input); %显示一下路径用的
%test_leachRW(s_mobility,s_input);
%test_DCAMM(s_mobility,s_input);
% test_SEP(s_mobility,s_input);
% test_DEEC(s_mobility,s_input);

test_throughput_DCAMM(s_mobility,s_input);

