function [value,isterminal,direction] = MyEvents(~,x)

g = max(x(1:end/3));
value = g-200;
isterminal = 1;
direction = 0;

end