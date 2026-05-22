function u = resample(u, n)

if ( isempty(u) )
    return
end

for k = 1:numel(u.components)
    u.components{k} = resample(u.components{k}, n);
end

end
