function [value, isterminal, direction] = events(~, y)
        value = y(2);
        isterminal = 1;
        direction = 0;
end

