%% 
% Modelling Quasi Birth-and-Death Process
% 
% $$\left\lbrack \begin{array}{ccccc}B_{0\text{ }}  & B_{1\text{ }}  &  &  
% & \\C & A & D &  & \\ & C & A & D & \\ &  & \ddots  & \ddots  & \ddots \end{array}\right\rbrack 
% \left\lbrack \begin{array}{c}{\mathbf{x}}_0 \\{\mathbf{x}}_1 \\{\mathbf{x}}_2 
% \\\vdots \end{array}\right\rbrack =\left\lbrack \begin{array}{c}0\\0\\0\\\vdots 
% \end{array}\right\rbrack$$
% 
% where 
% 
% $B_0 =-0\ldotp 2I_2 ,B_{1\text{ }} =0\ldotp 3I_{2\text{ }} ,C=0\ldotp 2I_2 
% \ldotp D=0\ldotp 3I_2 \text{ }\mathrm{and}\text{ }A=\left\lbrack \begin{array}{cc}-0\ldotp 
% 6 & 0\ldotp 2\\0\ldotp 1 & -0\ldotp 7\end{array}\right\rbrack$.
% 
% Assuming ${\mathbf{x}}_{\mathbf{i}+1} ={R\text{ }\mathbf{x}}_{\mathbf{i}}$. 
% Then, we have 
% 
% $\left\lbrace \begin{array}{cc}\left(B_0 +B_1 R\right){\mathbf{x}}_0 =0 
% & \left(1\right)\\\left(C+\mathit{AR}+{\mathit{DR}}^2 \right){\mathbf{x}}_0 
% =0 & \left(2\right)\end{array}\right.$.
% 
% *Defining variables:*

n = 2;
B0 = -0.2*eye(n)
B1 = 0.3*eye(n)
C = 0.2*eye(n)
D = 0.3*eye(n)
A = [-0.6, 0.2; 0.1, -0.7]
%% 
% *Iterative schemes:*
% 
% Since $\mathrm{det}\left(A\right)\not= 0$, 
% 
% $$\begin{array}{l}C+\text{AR}+{\text{DR}}^2 =0\\-\mathit{AR}=C+{\mathit{DR}}^2 
% \\R=-A^{-1} \left(C+{\mathit{DR}}^2 \right)\\R_{k+1} =-A^{-1} \left(C+{\mathit{DR}}_k^2 
% \right)\end{array}$$

mulA = -1/(-0.6*-0.7 - 0.1*0.2)*[-0.7, -0.2; -0.1, -0.6];
cons = mulA*C;
coef = mulA*D;
R = 0.333*eye(n);
maxIter = 100;
%Solve matrix quadratic eqaution
for i = 1:maxIter
    newR = cons + coef*R^2;
    if i == maxIter
        fprintf('Error matrix:\n')
        abs(newR - R)
        fprintf('Numerical output:\n')
        newR
    end 
    R = newR;
end 
%Solve for x
precision = 10;
coefMat = round(B0 + B1*R, precision)
rref(coefMat) % Numerically unstable without rounding the coefficient matrix
%% 
% Therefore, 
% 
% $${\mathbf{x}}_0 =\left\lbrack \begin{array}{c}2\\1\end{array}\right\rbrack 
% \ldotp$$
