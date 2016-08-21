function pass = fn_flares(pass, param)

% this functions counts the flares into the 'pass' structure
% summarizes data about their magnitudes, durations, etc.
% and saves them in 'pass.flares' field
%
% detailed description on https://habrahabr.ru/post/307212/

% initialization before the main loop
counter_flare = 0;
flare_is_now = 0;
clear flares
flares = struct('time',[],'dist',[],'alt',[],'azimuth',[],'inst_magn',[],'vis_magn',[],'dur',[]);

% main loop
for step = 1:length(pass.time)

    if pass.magnitude(step) < param.min_m
        % flare is present now
        if flare_is_now == 0

            % flare just began
            flare_is_now = 1;
            counter_flare = counter_flare + 1;
            clear current_flare
            current_flare.time_beg = pass.time(step);
            % integrate the energy of the flare, [J/m^2]
            if step < length(pass.time)
                current_flare.energy = pass.flux(step)*...
                    (pass.time(step+1) - pass.time(step));
            else
                % this was the last data point, nothing to add
            end

            % create field 'current_flare.brightest' 
            % with the data about the brightest point of the flare
            current_flare.brightest.magn = pass.magnitude(step);
            current_flare.brightest.time = pass.time(step);
            current_flare.brightest.azimuth = pass.azimuth(step);
            current_flare.brightest.alt = pass.alt(step);
            current_flare.brightest.distance = pass.distance(step);
        else
            % flare still goes on
            % increment the energy of the flare
            if step < length(pass.time)
                current_flare.energy = current_flare.energy + pass.flux(step)*...
                    (pass.time(step+1) - pass.time(step));
            else
                % this was the last data point, nothing to add
            end
            % if here flare is brighter, renew the data about the brightest point
            if pass.magnitude(step) < current_flare.brightest.magn
                current_flare.brightest.magn = pass.magnitude(step);
                current_flare.brightest.time = pass.time(step);
                current_flare.brightest.azimuth = pass.azimuth(step);
                current_flare.brightest.alt = pass.alt(step);
                current_flare.brightest.distance = pass.distance(step);
            end
        end
    else
        % no flare now
        if flare_is_now == 0
            
            % still no flare
        else
            % flare just ended, sum up the results
            flare_is_now = 0;
            current_flare.time_end = pass.time(step-1);
            current_flare.duration = current_flare.time_end - current_flare.time_beg;
            if current_flare.duration > param.frame_duration
                % if the flare is long, we can see it directly
                current_flare.vis_magn = current_flare.brightest.magn;
            else
                % if the flare is short, one should divide its energy
                % over the exposure time of an eye/camera to get the
                % apparent magnitude
                current_flare.vis_magn = ...
                    min(-2.5*log10(current_flare.energy/param.frame_duration/param.flux_m0), param.min_m);
            end

            % save the data in the 'flares' structure
            flares.time(counter_flare,1) = current_flare.brightest.time;
            flares.dist(counter_flare,1) = current_flare.brightest.distance;
            flares.azimuth(counter_flare,1) = current_flare.brightest.azimuth;
            flares.alt(counter_flare,1) = current_flare.brightest.alt;
            flares.inst_magn(counter_flare,1) = current_flare.brightest.magn;
            flares.vis_magn(counter_flare,1) = current_flare.vis_magn;
            flares.dur(counter_flare,1) = current_flare.duration;
        end
    end

end

% return the result
flares.num_of_flares = counter_flare;
pass.flares = flares;

end

