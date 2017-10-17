n = 3;
%angle
theta = [0, pi/4, pi/4];
%rate of change of the end effector coordinates 
%with respect to the joint coordinates
di = ADfun('direct_kinematics',n);
[f, A] = feval(di,theta);

%for given tip velocity x, choose a corresponding set of joint velocities
%is to choose the joint velocities that minimizes 2-norm
x = [0;1];
min_norm = A'*inv(A*A')*x;

%In many cases, it is more costly to move some joins than other.
%so, instead of minimizing 2-norm, we fid a solution which minimizes a
%different norm
D = diag([6,3,1]);
norm_d = inv(D)*A'*inv(A*inv(D)*A')*x;