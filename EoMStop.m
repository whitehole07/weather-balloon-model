function [value, isterminal, direction] = EoMStop(t, x)
value      = (t~=0 && (x(1) <= 0 && x(2) <= 0));
isterminal = 1;   % Stop the integration
direction  = 0;