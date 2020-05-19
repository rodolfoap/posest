%POSEST Estimate pose from noisy 3D-2D correspondences.
%
% [rt, idxOutliers, focal]=posest(pts2D, pts3D, inlPcent, K, NLrefine, verbose);
%        pts2D, pts3D are the matched 2D-3D point coordinates
%        inlPcent is the expected percentage of inliers (>=0.5)
%        K is either
%          - the 3x3 camera intrinsic calibration matrix
%          - a 1x2 vector with the image width and height. The camera focal length is estimated using P4Pf
%        inlPcent is the expected fraction of inliers in the input pairs
%        NLrefine controls non-linear refinement, can be one of: 'norefine', 'repr_err', 'repr_err_mlsl', 'objspc_err'
%        or empty (implying 'norefine', default)
%        verbose is optional
%
% Returns the 6x1 rotation - translation vector rt, the indices of outlying point pairs from pts2D, pts3D
%        and the camera focal length.
%
% Notes::
%   MEX file.
%

if ~exist('posest', 'file')
  error('you need to build the MEX version of posest, see mex');
end
