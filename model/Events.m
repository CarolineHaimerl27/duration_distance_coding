function [value,isterminal,direction] = Events(~,x)
% g = mean(x(1:50));
% value = g-0.85;

% changed this back to original event function that you sent me per mail
% since this one did not work

index = length(x)+(-x(1):0); %(x(1) number of neurons)
firing_rate = x(index);
g = max(firing_rate);
value = g-50;
isterminal = 1;
direction = 0;

end