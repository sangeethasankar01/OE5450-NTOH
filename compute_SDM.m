function [I,J] = COMPUTE_SDM(X_cp,Y_cp,X_b,Y_b,phi,L)


% PURPOSE
% - Compute the integral expression for constant strength source panels
% - Source panel strengths are constant, but can change from panel to panel
% - Geometric integral for panel-normal    : I(ij)
% - Geometric integral for panel-tangential: J(ij)

% INPUTS
% - X_cp  : X-coordinate of control points
% - Y_cp  : Y-coordinate of control points
% - X_b  : X-coordinate of boundary points
% - Y_b  : Y-coordinate of boundary points
% - phi : Angle between positive X-axis and interior of the panel
% - L : Length of panel
% 
% OUTPUTS
% - I   : Value of panel-normal integral (Eq. 6.266 in Anderson or Ref [1])
% - J   : Value of panel-tangential integral (Eq. 6.267 in Anderson or Ref [2])

% Number of panels
con_pt = length(X_cp);                                                        % Number of panels/control points

% Initialize arrays
I = zeros(con_pt,con_pt);                                                   % Initialize I integral matrix
J = zeros(con_pt,con_pt);                                                   % Initialize J integral matrix
% Compute integral
for i = 1:1:con_pt                                                         % Loop over i panels
    for j = 1:1:con_pt                                                      % Loop over j panels
        if (j ~= i)                                                         % If the i and j panels are not the same
            % Compute intermediate values
            A  = -(X_cp(i)-X_b(j))*cos(phi(j))-(Y_cp(i)-Y_b(j))*sin(phi(j));      % A term
            B  = (X_cp(i)-X_b(j))^2+(Y_cp(i)-Y_b(j))^2;                           % B term
            Cn = sin(phi(i)-phi(j));                                        % C term (normal)
            Dn = -(X_cp(i)-X_b(j))*sin(phi(i))+(Y_cp(i)-Y_b(j))*cos(phi(i));      % D term (normal)
            Ct = -cos(phi(i)-phi(j));                                       % C term (tangential)
            Dt = (X_cp(i)-X_b(j))*cos(phi(i))+(Y_cp(i)-Y_b(j))*sin(phi(i));       % D term (tangential)
            E  = sqrt(B-A^2);                                               % E term
            if (~isreal(E))
                E = 0;
            end
            
            % Compute I (needed for normal velocity), Ref [1]
            term1  = 0.5*Cn*log((L(j)^2+2*A*L(j)+B)/B);                     % First term in I equation
            term2  = ((Dn-A*Cn)/E)*(atan2((L(j)+A),E) - atan2(A,E));        % Second term in I equation
            I(i,j) = term1 + term2;                                         % Compute I integral
            
            % Compute J (needed for tangential velocity), Ref [2]
            term1  = 0.5*Ct*log((L(j)^2+2*A*L(j)+B)/B);                     % First term in J equation
            term2  = ((Dt-A*Ct)/E)*(atan2((L(j)+A),E) - atan2(A,E));        % Second term in J equation
            J(i,j) = term1 + term2;                                         % Compute J integral
        end
        
        % Zero out any NANs, INFs, or imaginary numbers
        if (isnan(I(i,j)) || isinf(I(i,j)) || ~isreal(I(i,j)))
            I(i,j) = 0;
        end
        if (isnan(J(i,j)) || isinf(J(i,j)) || ~isreal(J(i,j)))
            J(i,j) = 0;
        end
    end
end
