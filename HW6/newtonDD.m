% Newton Divided Difference Interpolation

function newtonDD(x, y)
    n = size(x,1);
    n2 = size(x,2);
    if n < n2
        x = x';
        y = y';
        n = n2;
    end
    
    coefs = newtonDDsetup(x, y);
    
    m = 10*n;
    minx = min(x);
    maxx = max(x);
    t = minx:(maxx-minx)/m:maxx;
    
    val = newtonEval(t, coefs, x);
    
    hold off;
    plot(t,val);
    hold on;
    for i=1:1:n
        plot(x(i),y(i),'k*');
    end
end
    
    
% Evaluate divided difference interpolant
function value = newtonEval(t, coefs, x)
    n = size(x,1);
    m = size(t,1);
    for i=1:1:m
        value(i) = coefs(n);
    end
    for i=n-1:-1:1
        value = value .* (t - x(i)) + coefs(i);
    end
end


% Set up divided difference coefficients
function coefs = newtonDDsetup(x, y)
    n = size(x,1);

    % DD level 0
    for i=1:1:n
        coefs(i) = y(i);
    end

    % DD higher levels (bottom to top, overwrite lower entries as they are finished)
    for level=2:1:n
        for i=n:-1:level
            dx = x(i) - x(i-level+1);
            coefs(i) = ( coefs(i)-coefs(i-1) ) / dx;
        end
    end
end
