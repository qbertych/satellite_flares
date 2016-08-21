function pass = fn_reflection (param, reflection_intervals)

    % this function should be called for a time interval 
    % when the satellite is above horizon and not in the shadow
    % since it does not verify it
    %
    % it calculates the coordinates as well as brightness
    % of the flares with a variable time step
    % and returns the structure named 'pass'
    %
    % detailed description on https://habrahabr.ru/post/307212/

    % rename used variables
    % world
    rad_earth = param.rad_earth;
    GM = param.GM;
    omega_obs = param.omega_obs;
    atmosphere_transmission = param.atmosphere_transmission;
    flux_sun = param.flux_sun;
    flux_m0 = param.flux_m0;
    min_m = param.min_m;
    k_sun = param.k_sun;
    % satellite orbit
    h_orbit = param.h_orbit;
    incl_sat = param.incl_sat;
    theta_asc = param.theta_asc;
    alpha_sat_0 = param.alpha_sat_0;
    % satellite mirror
    d_rot = param.d_rot;
    omega_rot = param.omega_rot;
    n_mirror_0 = param.n_mirror_0;
    gauss_w = param.gauss_w;
    sq_mirror = param.sq_mirror;
    % observer coordinates
    phi_obs = param.phi_obs;
    theta_obs_0 = param.theta_obs_0;
    % simulation parameters
    dt2 = param.dt2;
    dt3 = param.dt3;
    dt4 = param.dt4;
    
    
    % normalize all the vectors
    k_sun = fn_norm(k_sun);
    d_rot = fn_norm(d_rot);
    n_mirror_0 = fn_norm(n_mirror_0);

    % calculate some aux values
    rad_orbit = rad_earth + h_orbit;                    % radius of the orbit from the center of Earth
    r_np = [0 0 rad_earth];                             % coordinates of the North pole
    dir_an = [cosd(theta_asc), sind(theta_asc), 0];       % direction on the acsending node
    dir_b = [-sind(theta_asc)*cosd(incl_sat),...
        cosd(theta_asc)*cosd(incl_sat), sind(incl_sat)];   % direction on the apsis (highest latitude)
    omega_sat = sqrt(2*GM/rad_orbit)/(2*pi*rad_orbit);  % satellite angular velocity, s^-1
    angle_drot_nmirror = acosd(fn_dot(d_rot,n_mirror_0));  % angle between the mirror and the sat rotation axis
    n_mir_parallel = d_rot*fn_dot(d_rot,n_mirror_0);    % component of n_mirror_0 parallel to d_rot, not normalized!
    n_mir_perp0 = n_mirror_0 - n_mir_parallel;          % component of n_mirror_0 perpendicular to d_rot
    n_mir_perp1 = fn_cross(d_rot, n_mir_perp0);         % one more perpendicular to d_rot
    n_mir_perp1 = n_mir_perp1...
        / sqrt(fn_dot(n_mir_perp1,n_mir_perp1))...
        * sqrt(fn_dot(n_mir_perp0,n_mir_perp0));        % normalize to the same length
    flux_min = flux_m0*10^(-min_m/2.5);                 % minimal flux, below which reflection is considered to be invisible
    dir_obs_z = sind(phi_obs);                           % z-coordinate of the observer, does not change
    angle_div = 3*gauss_w;                              % angle limits for the reflection cone

    
    % initialization before the main loop
    step = 1;
    pass = struct('time',[],'distance',[],'alt',[],'azimuth',[],'flux',[],'magnitude',[]);
    
    for counter_intervals = 1:size(reflection_intervals,1)

        % beg and end of current interval when sat is not inshadow
        time_beg = reflection_intervals(counter_intervals,1);
        time_end = reflection_intervals(counter_intervals,2);
        current_time = time_beg;
        
        % calculate reflections within given interval
        while current_time < time_end

            % current time
            pass.time(step,1) = current_time;

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
            distance = sqrt(fn_dot(r_sat_obs, r_sat_obs));
            pass.distance(step,1) = distance;

            % basis on the skychart 
            r_obs_np = r_np - r_obs;                                        % vector from the observer to the North pole
            dist_obs_np = sqrt(sum(r_obs_np.^2));
            cos_angle_obs_np = fn_dot(r_obs_np, -dir_obs)/dist_obs_np;
            dir_nord = r_obs_np + cos_angle_obs_np*dist_obs_np*dir_obs;
            dir_nord = fn_norm(dir_nord);                                   % vector to the nord
            dir_west = fn_cross(dir_obs, dir_nord);
            dir_west = fn_norm(dir_west);                                   % vector to the west

            % coordinates of the satellite on the sky map
            pass.alt(step,1) = 90 - acosd(fn_dot(-dir_sat_obs, dir_obs));
            map_proj_nord = fn_dot(-dir_sat_obs, dir_nord);
            map_proj_west = fn_dot(-dir_sat_obs, dir_west);
            pass.azimuth(step,1) = mod(-atan2d(map_proj_west, map_proj_nord), 360);

            % calculate optimal direction of the mirror,i.e. mirror oriented 
            % like this will reflect light exactly on the observer
            n_opt = fn_norm((dir_sat_obs - k_sun)/2);

            % angle between d_rot and n_opt
            angle_drot_nopt = acosd(fn_dot(d_rot, n_opt));

            % CHECK 3: if angle (d_rot, n_opt) is close to angle (d_rot, n_mir0), proceed
            % this means that geometry can in principle give a reflection 
            if (angle_drot_nopt > (angle_drot_nmirror-angle_div))...
                    && (angle_drot_nopt < (angle_drot_nmirror+angle_div))
            
                % current orientation of the mirror
                n_mir = n_mir_parallel...
                    + n_mir_perp0*cosd(360*omega_rot*current_time)...
                    + n_mir_perp1*sind(360*omega_rot*current_time);

                % angle between current orientation of the mirror and n_opt
                angle_nmir_nopt = acosd(fn_dot(n_mir, n_opt));

                % CHECK 4: if n_mir is close to n_opt, proceed and calculate the reflection
                if angle_nmir_nopt < angle_div

                    % calculate the angle between the direction on the observer
                    % and direction of the main reflection
                    cos_angle_nmir_ksun = fn_dot(n_mir, -k_sun);
                    k_refl = 2*cos_angle_nmir_ksun*n_mir + k_sun;
                    cos_angle_refl_obs = abs(fn_dot(k_refl, dir_sat_obs));
                    angle_refl_obs = acosd(cos_angle_refl_obs);

                    % flux at the observer's place
                    refl_fraction = 2/(pi*degtorad(gauss_w)^2) * exp (-2*angle_refl_obs^2/gauss_w^2);     % now it is Gaussian
                    int_reflection = flux_sun * sq_mirror * cos_angle_nmir_ksun;
                    flux = max (int_reflection * refl_fraction ... 
                             * atmosphere_transmission / distance^2, flux_min);   % W/m^2

                    % save results
                    pass.flux(step,1) = flux;
                    pass.magnitude(step,1) = min(-2.5*log10(flux/flux_m0), min_m);
                    current_time = current_time + dt4;
                else
                    % observer is far from reflection, but might reach it soon
                    pass.flux(step,1) = flux_min;
                    pass.magnitude(step,1) = min_m;
                    current_time = current_time + dt3;
                end
            else
                % observer is far from reflection
                pass.flux(step,1) = flux_min;
                pass.magnitude(step,1) = min_m;
                current_time = current_time + dt2;
            end

            step = step+1;
        end

    end

end
