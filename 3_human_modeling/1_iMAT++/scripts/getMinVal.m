function [val_f,type_f, val_r, type_r] = getMinVal(OFD, epsilon_f,epsilon_r)

if OFD > 1e-5 % 1e-5 is the tol_flux, see supplement text for details
    val_f = -1;%no need for FVA ==> level 1
    type_f = 'NoNeed'; % no need
elseif OFD > 1e-7 % 1e-5 ~ 1e-7; test if it is a consistent small
    val_f = max(epsilon_f - 1e-5,1e-5); 
    type_f = 'smallFlux';% 1e-7 is the tol_zero, see supplement text for details
else
    val_f = max(epsilon_f - 1e-5,1e-5);
    type_f = 'normal';
end

if -OFD > 1e-5
    val_r = -1;%no need for FVA ==> level 1
    type_r = 'NoNeed'; % no need
elseif -OFD > 1e-7 % test if it is a consistent small
    val_r = 1e-7;
    type_r = 'smallFlux';
else
    val_r = max(epsilon_r - 1e-5,1e-5);
    type_r = 'normal';
end

end