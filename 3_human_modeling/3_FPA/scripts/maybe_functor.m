function c = maybe_functor(f, a, b)
    
    if isnan(a) && isnan(b)
        c = nan;
    elseif ~isnan(a) && isnan(b)
        c = a;
    elseif isnan(a) && ~isnan(b)
        c = b;
    else 
        c = f(a,b);
    end
end
