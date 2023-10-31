function test_Animate(s_mobility,s_input)
    v_t = 0:s_input.TIME_STEP:s_input.SIMULATION_TIME;   % 0 0.1 0.2 0.3 0.4...
    
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

    title(cat(2,'Simulation time (sec): ',num2str(s_mobility.SIMULATION_TIME)));
    xlabel('X (meters)');
    ylabel('Y (meters)');
    title('Radom Waypoint mobility');
    %ht在工作区是Text,标注
    ht = text(s_input.V_POSITION_X_INTERVAL(1),s_input.V_POSITION_Y_INTERVAL(2),cat(2,'Time (sec) = 0'));
    %ht = text(min(vs_node(1).v_x),max(vs_node(1).v_y),cat(2,'Time (sec) = 0'));
    %伸缩坐标轴
    axis([s_input.V_POSITION_X_INTERVAL(1) s_input.V_POSITION_X_INTERVAL(2) s_input.V_POSITION_Y_INTERVAL(1) s_input.V_POSITION_Y_INTERVAL(2)]);
    %axis([min(vs_node(1).v_x) max(vs_node(1).v_x) min(vs_node(1).v_y) max(vs_node(1).v_y)]);
    hold off;

    for timeIndex = 1:length(v_t)
        %每个时刻
        t = v_t(timeIndex);
        set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
        %%%%%%% LEACH begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test_leachRW(s_mobility,s_input,vh_node_pos,vs_node,timeIndex);
        %%%%%%% LEACH end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nodeIndex = 1:s_mobility.NB_NODES
            set(vh_node_pos(nodeIndex),'XData',vs_node(nodeIndex).v_x(timeIndex),'YData',vs_node(nodeIndex).v_y(timeIndex));
        end
        drawnow;
    end
end