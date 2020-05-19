matchesfile='../examples/mono/32D.txt';
calibfile='../examples/mono/K.txt';

% read in matches
[x, y, z, u, v]=textread(matchesfile, '%f%f%f%f%f', 'commentstyle', 'shell');
pts3D(:, :)=[x'; y'; z']';
pts2D(:, :)=[u'; v']';
%pt3D, pts2D

% read intrinsics
K=dlmread(calibfile, ' ', [0 0 2 2]);

% estimate pose

%[rt, idxOutliers]=posest(pts2D, pts3D, 0.8, K, 'objspc_err'); % LHM
%[rt, idxOutliers]=posest(pts2D, pts3D, 0.8, K, 'repr_err_mlsl'); % MLSL (slow!)
[rt, idxOutliers]=posest(pts2D, pts3D, 0.8, K, 'repr_err');
%[rt, idxOutliers]=posest(pts2D, pts3D, 0.8, K, 'repr_err', 1); % verbose
nbOutliers=max(size(idxOutliers));

format long g
disp('Estimated pose (rt)');
rt', nbOutliers

% estimate pose + focal length; image dimensions retrieved from K
[rtf, idxOutliers, foc]=posest(pts2D, pts3D, 0.8, [K(1, 3)*2, K(2, 3)*2], 'repr_err');

disp('Estimated pose (rt) and focal length');
rtf', nbOutliers, foc
