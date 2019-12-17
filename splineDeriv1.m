function d1 = splineDeriv1(pp, xx, x)

n = length(x);
y = zeros(1,length(xx));

for (i = 1:n-1)

    dp = x(i);

    slope = pp.coefs(i,:);
    der(i,1) = slope(1) * 3;
    der(i,2) = slope(2) * 2;
    der(i,3) = slope(3);

    f = @(x) der(i,1)*(x-dp)^2 + der(i,2)*(x-dp) +  der(i,3);

    for (j = 1:length(xx))
        if ((xx(j) >= x(i) && xx(j) <= x(i+1)) || (xx(j) < x(1) || xx(j) > x(end)))
            y(j) = f(xx(j));
        end
    end
end

d1 = y;
