function plot_orbit_SPICE(target, l_star, varargin)
    frame    = 'ECLIPJ2000'; % ECLIPJ2000
    abcorr   = 'NONE';
    observer = 'Sun';   
    
    if strcmp(target, 'Earth')
        Period_days = 365;
    elseif strcmp(target, 'Mars')
        Period_days = 687;
    elseif strcmp(target, 'Vesta')
        Period_days = 1325;
    elseif strcmp(target, 'Ceres')
        Period_days = 1682;
    else
        disp('Target name might be wrong??')
    end
    
    etmin = cspice_str2et('January 01, 2000 12:00 AM UTC');
    %  etmin = cspice_str2et('2024 Aug 11 12:00:00')
    etmax = etmin + 60 * 60 * 24 * Period_days;
    et = linspace(etmin,etmax,floor((etmax-etmin)/(24*60*60))) ;  %every day
    starg = cspice_spkezr( target, et, frame, abcorr, observer);

    plot(starg(1,:)/l_star, starg(2,:)/l_star, varargin{:})
end