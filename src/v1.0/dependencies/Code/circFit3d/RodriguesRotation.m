function rotatedPts = RodriguesRotation(providedPts, origNormal, newNormal)
% rotate provided points based on a staring and ending vector
origNormal = origNormal/norm(origNormal);
newNormal = newNormal/norm(newNormal);
numProvidedPts = size(providedPts, 1);
theta = acos(dot(origNormal, newNormal) );
rotationNormal = cross(origNormal, newNormal);
rotationNormal = rotationNormal / norm(rotationNormal );
if ~sum(isnan(rotationNormal)),
    rotatedPts = providedPts.*(cos(theta)*ones(numProvidedPts, 3)) + ...
        cross(ones(numProvidedPts,1)*rotationNormal, providedPts, 2) .* (sin(theta)*ones(numProvidedPts, 3)) + ...
        (ones(numProvidedPts,1)*rotationNormal) .* (dot(ones(numProvidedPts,1)*rotationNormal, providedPts, 2)*ones(1, 3)) .* ...
        ((1-cos(theta))*ones(numProvidedPts, 3)); % Rodrigues' formula
else
    rotatedPts = providedPts;
end


return