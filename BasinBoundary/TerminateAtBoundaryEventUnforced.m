function [value, isterminal, direction] = TerminateAtBoundaryEventUnforced(t, y, Mrhs)
value = round(Mrhs.bI(y(1:2)'),0) - 2; %zero when at boundary or origin, one otherwise
isterminal = 1;   % Stop the integration
direction  = 0;
end

