%% Define boundary conditions
% Fixed - Free
BC.U.location = [x == 0];
BC.U.val = [0];

BC.V.location = [x == 0];
BC.V.val = [0];

BC.W.location = [x == 0];
BC.W.val = [0];

BC.Theta1.location = [x == 0];
BC.Theta1.val = [0];

BC.Theta2.location = [x == 0];
BC.Theta2.val = [0];

BC.Theta3.location = [x == 0];
BC.Theta3.val = [0];

%% Fixed - Fixed
BC.U.location = [x == 0, x == 1];
BC.U.val = [0, 0];

BC.V.location = [x == 0, x == 1];
BC.V.val = [0, 0];

BC.W.location = [x == 0, x == 1];
BC.W.val = [0, 0];

BC.Theta1.location = [x == 0, x == 1];
BC.Theta1.val = [0, 0];

BC.Theta2.location = [x == 0, x == 1];
BC.Theta2.val = [0, 0];

BC.Theta3.location = [x == 0, x == 1];
BC.Theta3.val = [0, 0];

%% 
% Fixed - Roller
BC.U.location = [x == 0];
BC.U.val = [0];

BC.V.location = [x == 0, x == 1];
BC.V.val = [0, 0];

BC.W.location = [x == 0, x == 1];
BC.W.val = [0, 0];

BC.Theta1.location = [x == 0, x == 1];
BC.Theta1.val = [0, 0];

BC.Theta2.location = [x == 0, x == 1];
BC.Theta2.val = [0, 0];

BC.Theta3.location = [x == 0];
BC.Theta3.val = [0];

%% 
% Pinned - Roller
BC.U.location = [x == 0];
BC.U.val = [0];

BC.V.location = [x == 0, x == 1];
BC.V.val = [0, 0];

BC.W.location = [x == 0, x == 1];
BC.W.val = [0, 0];

BC.Theta1.location = [x == 0, x == 1];
BC.Theta1.val = [0, 0];

BC.Theta2.location = [x == 0, x == 1];
BC.Theta2.val = [0, 0];

BC.Theta3.location = [NaN];
BC.Theta3.val = [NaN];