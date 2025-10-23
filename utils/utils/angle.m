% Compute the angle between two vectors.
function out = angle(vec1, vec2)
    out = acos(dot(vec1, vec2) / (norm(vec1) * norm(vec2)));
end
