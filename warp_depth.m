%   TGV2-L2 Depth Image Upsampling
%
%   Author: David Ferstl
%
%   If you use this file or package for your work, please refer to the
%   following papers:
% 
%   [1] David Ferstl, Christian Reinbacher, Rene Ranftl, Matthias Rüther 
%       and Horts Bischof, Image Guided Depth Upsampling using Anisotropic
%       Total Generalized Variation, ICCV 2013.
%
%   License:
%     Copyright (C) 2013 Institute for Computer Graphics and Vision,
%                      Graz University of Technology
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see
%     <http://www.gnu.org/licenses/>.
%
% WARP_DEPTH Warping of an depth image + intensity image into another
%   camera coordinate system
%   [ DTOF_WARPED, IR_WARPED ] = WARP_DEPTH( P_TOF, K_GRAY, DTOF, IR, M, N, FLAGS_TOF )
%
%   [IN]
%       P_TOF  ...    Projection Matrix of the depth camera 
%       K_GRAY ...    Camera Matrix of the second target camera
%       DTOF   ...    depth image
%       IR     ...    Intensity image
%       M,N    ...    Targe image size
%       FLAGS_TOF ... Optional parameter to set invalid depth pixels
%   [OUT]
%       DTOF_WARPED ... Warped depth image
%       IR_WARPED ... Warped intensity image
%

function [ dtof_warped, ir_warped ] = warp_depth( P_TOF, K_GRAY, dtof, ir, M, N, flags_tof )

    if(nargin == 6)
        flags_tof = ones(size(dtof));
    end

    [Mt, Nt] = size(dtof);

    % Get ToF camera center out of projection matrix
    C_TOF = null(P_TOF);
    C_TOF = C_TOF./C_TOF(4);

    [xxt, yyt] = meshgrid(1:Nt, 1:Mt);
    ptst = [xxt(:)'; yyt(:)'; ones(1,Mt*Nt); zeros(1,Mt*Nt)];
    idxtof = sub2ind([Mt, Nt], ptst(2,:), ptst(1,:));

    % mapping of depth image points into 3D space
    lines3D = inv([P_TOF;0,0,0,1])*ptst;
    nor = 1 ./ sqrt(sum(lines3D.^2, 1));
    lines3D = bsxfun(@rdivide, lines3D, nor);
    ptsTof_real = repmat(C_TOF(1:3,:), [1, size(lines3D,2)]) + lines3D(1:3,:).*repmat(dtof(:)', [3,1]);
    ptsTof_real = ptsTof_real(:,flags_tof(:)' == 1);
    val_idxtof = idxtof(flags_tof(:)' == 1);

    % remap points into image space of intensity camera
    pts2D = K_GRAY*ptsTof_real;
    pts2D = pts2D ./ repmat(pts2D(3,:), [3,1]);

    d = ptsTof_real(3,:);
    img = zeros(M,N); img_idxtof = img;
    [~, sortD] = sort(d, 'descend');

    % fill target image with remapped points
    for i = 1:size(sortD,2)
        pt = round(pts2D(1:2,sortD(i)));
        if(pt(2) <= M && pt(2) > 0 && pt(1) <= N && pt(1) > 0)
            img(pt(2), pt(1)) = d(sortD(i));
            img_idxtof(pt(2), pt(1)) = ir(val_idxtof(sortD(i)));
        end
    end

    dtof_warped = img(1:M, 1:N);
    ir_warped = img_idxtof(1:M, 1:N);

end

