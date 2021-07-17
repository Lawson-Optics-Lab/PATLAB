function c = speedSoundWater(T)
%SPEEDSOUNDWATER Calculate the speed of sound in water with temperature.
%
% DESCRIPTION:
%     This function will be deprecated in a future version of k-Wave. Use
%     waterSoundSpeed instead.
%
% ABOUT:
%     author        - Bradley E. Treeby
%     date          - 11th August 2008
%     last update   - 4th April 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2008-2019 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% display warning
disp('WARNING: speedSoundWater will be deprecated in a future version of k-Wave. Please changes usage to waterSoundSpeed.')

% pass inputs to waterSoundSpeed
c = waterSoundSpeed(T);