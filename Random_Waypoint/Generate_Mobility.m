function s_mobility = Generate_Mobility(s_input)
    %The Random Waypoint mobility model.
    global s_mobility_tmp;
    global nodeIndex_tmp;
    
    s_mobility.NB_NODES = s_input.NB_NODES;
    s_mobility.SIMULATION_TIME = s_input.SIMULATION_TIME;
    s_mobility.OBSTACLE_EDGE = s_input.OBSTACLE_EDGE;

    %产生障碍区域
    for obstacle = 1:s_input.OBSTAClE_NUM
        s_mobility.OBSTACLE(obstacle).V_POSITION_X = randi([s_input.OBSTACLE_EDGE /2, s_input.V_POSITION_X_INTERVAL(2) - s_input.OBSTACLE_EDGE /2]);
        s_mobility.OBSTACLE(obstacle).V_POSITION_Y = randi([s_input.OBSTACLE_EDGE /2, s_input.V_POSITION_Y_INTERVAL(2) - s_input.OBSTACLE_EDGE /2]);
    end

    for nodeIndex_tmp = 1:s_mobility.NB_NODES  %每个节点循环一次
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Initialize:
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION = [];
        % 生成随机数
        task_num = s_input.TASK_NUM;
        randomNumber = ceil(rand() * task_num);
        % 调整范围
        if randomNumber > task_num
            randomNumber = task_num;
        end
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).TASK_ID = randomNumber;

        previousX = unifrnd(s_input.V_POSITION_X_INTERVAL(1),s_input.V_POSITION_X_INTERVAL(2));  %unifrnd开始random，这里是random初始位置
        previousY = unifrnd(s_input.V_POSITION_Y_INTERVAL(1),s_input.V_POSITION_Y_INTERVAL(2));
        previousDuration = 0;
        previousTime = 0;
        Out_setRestrictedWalk_random_waypoint(previousX,previousY,previousDuration,previousTime,s_input); %限制

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Promenade     
        while (s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end) < s_input.SIMULATION_TIME)
            if (s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING(end) == false)%现在是重新移动的时候了
                previousX = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end);
                previousY = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end);
                previousDuration = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end);
                previousTime = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end);
                Out_setRestrictedWalk_random_waypoint(previousX,previousY,previousDuration,previousTime,s_input);
            else
                %%%%%%%%Node is taking a pause:
                previousDirection = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(end);
                previousSpeed = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(end);
                previousX = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end);
                previousY = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end);
                previousTime = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end);
                previousDuration = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end);
                distance = previousDuration*previousSpeed;
                %%%
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end+1,1) = previousTime + previousDuration;
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end+1,1) = (previousX + distance*cosd(previousDirection));
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end+1,1) = (previousY + distance*sind(previousDirection));
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(end+1,1) = 0;
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(end+1,1) = 0;
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING(end+1,1) = false;
                s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end+1,1) = Out_adjustDuration_random_waypoint(s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end),unifrnd(s_input.V_PAUSE_INTERVAL(1),s_input.V_PAUSE_INTERVAL(2)),s_input);
            end
        end
        %%%%%%%%%%%%%%%%%%To have speed vectors as well rather than
        %%%%%%%%%%%%%%%%%%only the scalar value:
        nb_speed = length(s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE);
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_X = zeros(nb_speed,1);
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_Y = zeros(nb_speed,1);
        for s = 1:nb_speed
            speed = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(s);
            direction = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(s);
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_X(s) = speed*cosd(direction);
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_Y(s) = speed*sind(direction);
        end

        %%%%%%%%%%%%%%%%%%To remove null pauses:
        v_index = s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(1:end-1) == 0;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_X(v_index) = [];
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_Y(v_index) = [];

        %%%%%%%%%%%%%%%%%%To remove the too small difference at the end, if
        %%%%%%%%%%%%%%%%%%there is one:
        if ((s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end) - s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end-1)) < 1e-14)
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_X(end) = [];
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_Y(end) = [];
        end
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end) = s_input.SIMULATION_TIME;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end) = 0;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(end) = 0;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_X(end) = 0;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_Y(end) = 0;
        
        s_mobility.VS_NODE(nodeIndex_tmp) = struct('V_TIME',s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME,...
                                              'V_POSITION_X',s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X,...
                                              'V_POSITION_Y',s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y,...
                                              'V_SPEED_X',s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_X,...
                                              'V_SPEED_Y',s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_Y,...
                                              'TASK_ID',s_mobility_tmp.VS_NODE(nodeIndex_tmp).TASK_ID);

    end
    
    clear s_mobility_tmp;
    clear nodeIndex_tmp;
end

function Out_setRestrictedWalk_random_waypoint(previousX,previousY,previousDuration,previousTime,s_input)

    global s_mobility_tmp;
    global nodeIndex_tmp;

    x_tmp = previousX;
    y_tmp = previousY;
    time_tmp = previousTime + previousDuration;
    duration_tmp = Out_adjustDuration_random_waypoint(time_tmp,unifrnd(s_input.V_WALK_INTERVAL(1),s_input.V_WALK_INTERVAL(2)),s_input); %这里是random walktime
    direction_tmp = unifrnd(s_input.V_DIRECTION_INTERVAL(1),s_input.V_DIRECTION_INTERVAL(2)); %这里是random direction
    speed = unifrnd(s_input.V_SPEED_INTERVAL(1),s_input.V_SPEED_INTERVAL(2));                 %这里是random speed
    distance_tmp = speed*duration_tmp;
    if (distance_tmp == 0)
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end+1,1) = time_tmp;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end+1,1) =  x_tmp;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end+1,1) =  y_tmp;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(end+1,1) = direction_tmp;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(end+1,1) = speed;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING(end+1,1) = true;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end+1,1) = duration_tmp;
    else
        %The loop begins:
        flag_mobility_finished = false;
        while (~flag_mobility_finished)
            x_dest = x_tmp + distance_tmp*cosd(direction_tmp);
            y_dest = y_tmp + distance_tmp*sind(direction_tmp);
            flag_mobility_was_outside = false;
            if (x_dest > s_input.V_POSITION_X_INTERVAL(2))
                flag_mobility_was_outside = true;
                new_direction = 180 - direction_tmp;
                x_dest = s_input.V_POSITION_X_INTERVAL(2);
                y_dest = y_tmp + diff([x_tmp x_dest])*tand(direction_tmp);  
            end
            if (x_dest < s_input.V_POSITION_X_INTERVAL(1))
                flag_mobility_was_outside = true;
                new_direction = 180 - direction_tmp;
                x_dest = s_input.V_POSITION_X_INTERVAL(1);
                y_dest = y_tmp + diff([x_tmp x_dest])*tand(direction_tmp);
            end
            if (y_dest > s_input.V_POSITION_Y_INTERVAL(2))
                flag_mobility_was_outside = true;
                new_direction = -direction_tmp;
                y_dest = s_input.V_POSITION_Y_INTERVAL(2);
                x_dest = x_tmp + diff([y_tmp y_dest])/tand(direction_tmp); 
            end
            if (y_dest < s_input.V_POSITION_Y_INTERVAL(1))
                flag_mobility_was_outside = true;
                new_direction = -direction_tmp;
                y_dest = s_input.V_POSITION_Y_INTERVAL(1);
                x_dest = x_tmp + diff([y_tmp y_dest])/tand(direction_tmp);
            end
            current_distance = abs(diff([x_tmp x_dest]) + 1i*diff([y_tmp y_dest]));
            current_duration = Out_adjustDuration_random_waypoint(time_tmp,current_distance/speed,s_input);
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end+1,1) = time_tmp;
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end+1,1) = x_tmp;
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end+1,1) = y_tmp;
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DIRECTION(end+1,1) = direction_tmp;
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_SPEED_MAGNITUDE(end+1,1) = speed;
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_IS_MOVING(end+1,1) = true;
            s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_DURATION(end+1,1) = current_duration;
            if(flag_mobility_was_outside)
                time_tmp = time_tmp + current_duration;
                duration_tmp = duration_tmp - current_duration;
                distance_tmp = distance_tmp - current_distance;
                x_tmp = x_dest;
                y_tmp = y_dest;
                direction_tmp = new_direction;
            else
                flag_mobility_finished = true;
            end
        end
        %the loop ended
    end
end

%
function duration = Out_adjustDuration_random_waypoint(time,duration,s_input)
    if ((time+duration) >= s_input.SIMULATION_TIME)
        duration = s_input.SIMULATION_TIME - time;
    end
end