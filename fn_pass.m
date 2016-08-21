function pass_path = fn_pass(param)

% calculates if the next turn around the Earth will be above observer's
% horizon and not in the shadow
%
% if the turn is above horizon, the function returns the sky chart coordinates of the pass
% if this pass is not in the shadow, the function returns time interval when it is visible
% if satellite is still visible at the end of the turn, 
% the function works till the satellite will go behind horizon,
% but not more than for one more turn
%
% detailed description on https://habrahabr.ru/post/307212/

    % rename used variables   
    % world
    rad_earth = param.rad_earth;        
    GM = param.GM;                      
    omega_obs = param.omega_obs;        
    k_sun = param.k_sun;                
    % satellite orbit
    h_orbit = param.h_orbit;            
    incl_sat = param.incl_sat;          
    theta_asc = param.theta_asc;        
    alpha_sat_0 = param.alpha_sat_0;    
    % observer coordinates
    phi_obs = param.phi_obs;            
    theta_obs_0 = param.theta_obs_0;    
    % simulation parameters
    dt1 = param.dt1;                    

    
    % normalize all the vectors
    k_sun = fn_norm(k_sun);
    
    % calculate some aux values
    rad_orbit = rad_earth + h_orbit;                    % radius of the orbit from the center of Earth
    r_np = [0 0 rad_earth];                             % coordinates of the North pole
    dir_an = [cosd(theta_asc), sind(theta_asc), 0];       % direction on the acsending node
    dir_b = [-sind(theta_asc)*cosd(incl_sat),...
        cosd(theta_asc)*cosd(incl_sat), sind(incl_sat)];   % direction on the highest point of the sat
    omega_sat = sqrt(2*GM/rad_orbit)/(2*pi*rad_orbit);  % satellite angular velocity on the orbit, s^-1
    dir_obs_z = sind(phi_obs);                           % z-coordinate of the observer, does not change
    
    % initialization before the main loop
    current_time = 0;
    time_max = 1/omega_sat;
    flag_currently_passing = 0;                          % flag showing that the sat is above horizon
    flag_pass_ended = 0;                                 % flag showing that the sat just disappeared beyond the horizon
    pass_path = struct('time',[],'distance',[],'alt',[],'azimuth',[],'in_shadow',[]);
    
    % STEP 1: calculate if the sat is above horizon and in shadow
    while ((current_time < time_max) || (flag_currently_passing == 1)) && (current_time < 2*time_max) && (flag_pass_ended == 0)

        % coordinates of the observer
        theta_obs = mod(theta_obs_0 + 360*omega_obs*current_time, 360);
        dir_obs = [cosd(phi_obs)*cosd(theta_obs), cosd(phi_obs)*sind(theta_obs), dir_obs_z];
        r_obs = rad_earth * dir_obs;

        % coordinates of the satellite
        theta_sat = mod(alpha_sat_0 + 360*omega_sat*current_time, 360);
        dir_sat = cosd(theta_sat)*dir_an + sind(theta_sat)*dir_b;
        r_sat = rad_orbit*dir_sat;
        
        % vector from the satellite to the observer
        r_sat_obs = r_obs-r_sat;                
        dir_sat_obs = fn_norm(r_sat_obs);  

        % CHECK1: check if the satellite is above the horizon
        if fn_dot(dir_obs, dir_sat_obs) < 0

            if flag_currently_passing == 0
                % satellite just raised above the horizon,
                flag_currently_passing = 1;
                pass_step = 1;
                sun_above_hor_beg = acosd(fn_dot(k_sun, dir_obs))-90;
            else
                % satellite was above the horizon before
                pass_step = pass_step + 1;
            end
            
            % save current time and distance from the observer to the satellite
            pass_path.time(pass_step,1) = current_time;
            pass_path.distance(pass_step,1) = sqrt(fn_dot(r_sat_obs, r_sat_obs));

            % basis on the skychart 
            r_obs_np = r_np - r_obs;                                    % vector from the observer to the North pole
            dist_obs_np = sqrt(sum(r_obs_np.^2));
            cos_angle_obs_np = fn_dot(r_obs_np, -dir_obs)/dist_obs_np;
            dir_nord = r_obs_np + cos_angle_obs_np*dist_obs_np*dir_obs;
            dir_nord = fn_norm(dir_nord);                               % vector to the nord
            dir_west = fn_cross(dir_obs, dir_nord);
            dir_west = fn_norm(dir_west);                               % vector to the west

            % coordinates of the satellite on the sky chart
            pass_path.alt(pass_step,1) = 90 - acosd(fn_dot(-dir_sat_obs, dir_obs));
            map_proj_nord = fn_dot(-dir_sat_obs, dir_nord);
            map_proj_west = fn_dot(-dir_sat_obs, dir_west);
            pass_path.azimuth(pass_step,1) = mod(-atan2d(map_proj_west, map_proj_nord), 360);


            % CHECK 2: check that the satellite is not in the shadow
            cos_angle_sat_ksun = fn_dot(dir_sat,k_sun);                  
            dist_sat_ksun = rad_orbit*sqrt(1-cos_angle_sat_ksun^2);     % distance from the sat to the central axis of shadow
            if (dist_sat_ksun > rad_earth)
                % satellite is not in the shadow
                pass_path.in_shadow(pass_step,1) = 0;  
            else
                if  fn_dot(dir_sat_obs, k_sun) < 0
                    % satellite is closer to the Sun than Earth
                    pass_path.in_shadow(pass_step,1) = 0; 
                else
                    % satellite is in the shadow
                    pass_path.in_shadow(pass_step,1) = 1;
                end
            end
            
        else
            % satellite is below the horizon        
            % if pass just ended, exit from the function
            if flag_currently_passing == 1
                flag_currently_passing = 0;
                flag_pass_ended = 1;
                sun_above_hor_end = acosd(fn_dot(k_sun, dir_obs))-90;
                pass_path.sat_above_hor = max(pass_path.alt);
            end
 
        end
        current_time = current_time + dt1;
    end
    
    
    % STEP 2: if the pass was above horizon,
    % select the regions when sat is above horizon and not in shadow
    pass_path.reflection_intervals = [];
    if isempty(pass_path.time) == 0
        
        flag_in_shadow = 1;
        num_regions_not_in_shadow = 0;
        for counter_step = 1:size(pass_path.time,1)
            
            if pass_path.in_shadow(counter_step) == 0
                % satellite is not in shadow
                if flag_in_shadow == 0  
                    % sat is still not in shadow
                else 
                    % sat was in shadow and just went out of it
                    flag_in_shadow = 0;
                    num_regions_not_in_shadow = num_regions_not_in_shadow + 1;
                    pass_path.reflection_intervals(num_regions_not_in_shadow,1) = pass_path.time(counter_step);
                end
            else
                % satellite is in shadow
                if flag_in_shadow == 0
                    % sat was not in shadow and just went into it
                    flag_in_shadow = 1;
                    pass_path.reflection_intervals(num_regions_not_in_shadow,2) = pass_path.time(counter_step);
                else
                    % sat is still in shadow
                end                
            end
            
        end
        
        % if at the end sat was still not in shadow,
        % time of the end of the last interval will be time of the last
        % point above the horizon
        if flag_in_shadow == 0
            pass_path.reflection_intervals(num_regions_not_in_shadow,2) = pass_path.time(counter_step);
        end  
    
    end
    
    % save data about the sun above horizon
    if exist('sun_above_hor_beg', 'var') && exist('sun_above_hor_end', 'var')
        pass_path.sun_above_hor = (sun_above_hor_beg+sun_above_hor_end)/2;
    end
    
    % if the satellite was visible for too long, send a message
    if current_time >= 2*time_max
        disp ('Satellite is visible for more than 1 turn');
    end
    
    
end
