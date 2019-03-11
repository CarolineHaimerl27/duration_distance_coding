function [value,isterminal,direction] = Events(~,x)

index = length(x)+(-x(1):0); %(x(1) number of neurons)
firing_rate = x(index);
g = max(firing_rate);
value = g-50;
isterminal = 1;
direction = 0;

end