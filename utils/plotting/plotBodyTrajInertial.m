%% plotBodyTrajInertial: function description
function [ah, RVs_body_I_out] = plotBodyTrajInertial(bodyname, et_range, primaryBody, fh, fprop, SF_l, fidelity)

if isempty(fh)
	fig = gcf;
else
	fig = figure(fh);
end
ax = fig.CurrentAxes;

if ~exist('fidelity', 'var')
	fidelity = 'spice';
end


switch fidelity
case 'spice'
	RVs_body_I_out = cspice_spkezr(bodyname, et_range, 'ECLIPJ2000', 'none', primaryBody);

    % calculate orbital elements
	GM 		= cspice_bodvrd('sun', 'GM', 1);
	COEvec0 = rv2COEvec(RVs_body_I_out(1:3,1), RVs_body_I_out(4:6,1), GM);
	n0 		= sqrt(GM/COEvec0(1)^3);

	RVs_body_I_plt = RVs_body_I_out;
	Tperiod = 2*pi/n0;
	t_his = et_range - et_range(1);
	RVs_body_I_plt(:,t_his > Tperiod) = [];

case 'Keplerian'
	RV0_body_I = cspice_spkezr(bodyname, et_range(1), 'ECLIPJ2000', 'none', primaryBody);

	% calculate orbital elements
	GM 		= cspice_bodvrd('sun', 'GM', 1);
	COEvec0 = rv2COEvec(RV0_body_I(1:3), RV0_body_I(4:6), GM);
	MA0  	= keplerInv(COEvec0(end), COEvec0(2), 1E-08, 'TA');
	n0 		= sqrt(GM/COEvec0(1)^3);
	RVs_body_I_out = getBodyOrbitState_kep(COEvec0, MA0, n0, et_range - et_range(1), @(COE) COE2rvvec(COE, GM));

	RVs_body_I_plt = RVs_body_I_out;
	Tperiod = 2*pi/n0;
	t_his = et_range - et_range(1);
	RVs_body_I_plt(:,t_his > Tperiod) = [];

end

if ~isfield(fprop, 'Flag_3D') || fprop.Flag_3D
	ah = plot3(RVs_body_I_plt(1,:)/SF_l,RVs_body_I_plt(2,:)/SF_l,RVs_body_I_plt(3,:)/SF_l,'k--');
else
	ah = plot(RVs_body_I_plt(1,:)/SF_l,RVs_body_I_plt(2,:)/SF_l,'k--');
end
end