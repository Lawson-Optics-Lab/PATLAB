function [centerLoc, circleNormal, radius] = CircFit3D(circLocs)
% extends circFit() to provide the best fit for a circle in R3.
% circLocs nx3: x, y, z coordinates

meanLoc = mean(circLocs);
numCurPts = length(circLocs);
movedToOrigin = circLocs - ones(numCurPts, 1)*meanLoc;
% plot3(curLocs(:, 2), curLocs (:, 3), curLocs (:, 4)); grid on;
% xlabel('x'); ylabel('y'); zlabel('z');
[U, s, V]= svd(movedToOrigin);
circleNormal = V(:, 3);
circleLocsXY = RodriguesRotation(movedToOrigin, circleNormal, [0, 0, 1]);
[xc, yc, radius] = circFit(circleLocsXY(:, 1), circleLocsXY(:, 2));
% centerLoc = [xc, yc, 0] + meanLoc;
centerLoc = RodriguesRotation([xc, yc, 0], [0, 0, 1], circleNormal) + meanLoc;
end