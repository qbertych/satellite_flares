%%%%%%%%%%%%%%%%%%     SIMULATION PARAMETERS     %%%%%%%%%%%%%%%%%%%%%%%%%%

% world
param.rad_earth = 6400 * 1E3;           % Earth radius, m
param.GM = 6.67*1E-11 * 5.972*1E24;     % grav. constant times Earth mass, m^3/s^2
param.omega_obs = 1/(60*60*24);         % Earth angular velocity, turns/s
param.k_sun = [-0.99, 0, 0.16];         % sunlight wavevector (22 Jun: [-0.92,0,-0.4]; 21 Dec: [-0.92,0,0.4], 15 Oct: [-0.99,0,0.16])
param.atmosphere_transmission = 0.7;    % atmosphere transmission
param.flux_sun = 1367;                  % Sun irradiant flux on the Earth orbit, W/m^2
param.flux_m0 = 2.5*1E-8;               % flux corresponding to the apparent magnitude = 0, W/m^2
param.min_m = 6;                        % minimal visible apparent magnitude

% observer
param.phi_obs = 55.75;                  % latitude, deg (Moscow: 55.75)
param.theta_obs_0 = 170;                % longitude, deg (noon: 0, midnight: 180)

% satellite orbit
param.incl_sat = 97.75;                 % inclination, deg (ISS: 51.6, SSO for 600 km: 97.75, one of Iridium: 88.6)
param.theta_asc = -7.5;                 % longitude of the ascending node, deg
param.alpha_sat_0 = 0;                  % angle between AN and the sat in the plane of orbit at t=0, deg
param.h_orbit = 600 * 1E3;              % orbit height, m

% satellite mirror
param.omega_rot = 1;                    % rotation speed of the satellite, turns/s
param.d_rot = [0 1 1];                  % rotation axis of the satellite
param.n_mirror_0 = [1 0 0];             % normal to the mirror at t=0
param.gauss_w = 17.2;                   % divergence angle of a normal distrubution, deg (Iridium: 0.3, Mayak: 17.2 (FWHM 20 deg))
param.sq_mirror = 3.8;                  % mirror square, m^2 (Iridium: 1.6, Mayak: 3.8)

% simulation parameters
param.dt1 = 1;                          % time step if the satellite is below horizon or in the shadow, s
param.dt2 = 1;                          % same if satellite is above horizon, but reflection is far away
param.dt3 = 1 * 1E-2;                   % same if satellite is above horizon, reflection is close
param.dt4 = 1 * 1E-3;                   % same if reflection is seen
param.frame_duration = 0.03;            % exposition time of a camera/eye, s

% data saving options
save_data = 0;
folder_name = 'd:\Mayak_sat\results';


         

%%%%%%%%%%%%%%%%%%%%%    HERE STARTS THE CODE    %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% detailed description on https://habrahabr.ru/post/307212/

% calculate trajectory of the pass
pass_path = fn_pass(param);

% proceed if the satellite is above the horizon and not always in shadow
if isempty(pass_path.reflection_intervals) == 0

    % calculate the reflections and magnitudes of flares
    pass = fn_reflection(param, pass_path.reflection_intervals);
    pass = fn_flares(pass, param);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%      IMAGES       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(pass_path.time) == 0
    
    % satellite passes above the horizon, plot the pass on the skychart
    figure(1)
    % r = 1 on the horizon
    polar(degtorad(pass_path.azimuth), 1-pass_path.alt/90, '-b')
    hold on
    polar(degtorad(pass_path.azimuth(logical(pass_path.in_shadow))),...
        1-pass_path.alt(logical(pass_path.in_shadow))/90, '.-k')
    hold off; view(-90, 90)

    if isempty(pass_path.reflection_intervals) == 0

        % sat is not in shadow, show centers of flares on the skychart
        figure(1)
        hold on
        % r = 1 on the horizon
        polar(degtorad(pass.flares.azimuth), 1-pass.flares.alt/90,'. r'); 
        hold off

        % dynamics
        figure(2)
        subplot(211); cla
        plot(pass.time, pass.magnitude, 'r', 'LineWidth', 1); hold on
        plot(pass.flares.time, pass.flares.vis_magn, '. b'); hold off
        set(gca, 'YLim',[-7, param.min_m], 'YDir','reverse'); ylabel ('Inst. magn.'); grid on
        title(['Max brightness ', num2str(min(pass.flares.vis_magn), '%0.1f')])
        subplot(212); cla
        plot(pass.flares.time, pass.flares.dur, '. b'); hold on
        plot([0 max(pass.time)],[param.frame_duration param.frame_duration], '--k', 'LineWidth', 1.5); hold off
        ylabel ('Flare duration, s'); grid on
        xlabel ('Pass time, s'); grid on
        set(gcf, 'PaperPosition', [0 0 9 15])
        set (gca, 'XLim', [0 max(pass.time)])
        linkaxes([subplot(211), subplot(212)], 'x');
    else
        % satellite is always in shadow
        figure(2); subplot(211); cla; title(''); subplot(212); cla;
        linkaxes([subplot(211), subplot(212)], 'x');
    end    
else
    % satellite is below the horizon
    figure(1); cla;
    figure(2); subplot(211); cla; title(''); subplot(212); cla;
    linkaxes([subplot(211), subplot(212)], 'x');
end



%%%%%%%%%%%%%%%%%%     SAVE FIGURES AND DATA      %%%%%%%%%%%%%%%%%%%%%%%%%

if save_data == 1
    
    if exist(folder_name, 'dir') == 0
        mkdir(folder_name);
    end
    file_basename = 'flare';
    saveas(figure(1), fullfile(folder_name, [file_basename, '_map.fig']));
    saveas(figure(2), fullfile(folder_name, [file_basename, '_dynamics.fig']));
    save(fullfile(folder_name, [file_basename, '_pass_path.mat']), 'pass_path');
    save(fullfile(folder_name, [file_basename, '_pass.mat']), 'pass');
    
end
