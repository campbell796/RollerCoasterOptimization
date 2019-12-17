function d2 = splineDeriv2(pp, xx, x)

n = length(x);
y = zeros(1,length(xx));

for (i = 1:n-1)
    
    dp = x(i);
    
    slope = pp.coefs(i,:);
    der(i,1) = slope(1) * 6;
    der(i,2) = slope(2) * 2;
    
    f = @(x) der(i,1)*(x-dp) + der(i,2);
    
    for (j = 1:length(xx))
        if ((xx(j) >= x(i) && xx(j) <= x(i+1)) || (xx(j) < x(1) || xx(j) > x(end)))
            y(j) = f(xx(j));
        end
    end
end

d2 = y;
