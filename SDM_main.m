clear;
clc;

%% Define physical variable
Uo = 1;                 % Freestream velocity
AoA  = 0;               % Angle of attack [deg]
N_Bp = 25;              % Number of boundary points (N+1)
tO   = (360/(N_Bp-1))/2;% Boundary point angle offset [deg]

% Plotting flags
flagPlot = [1;          % Shape polygon with panel normal vectors
            1;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Analytical and SPM pressure coefficient plot
            1;          % Streamlines
            1];         % Pressure coefficient contours

%% CREATE CIRCLE BOUNDARY POINTS

% Create angles to compute boundary points
theta = linspace(0,360,N_Bp)'; % Create angles for computing boundary point locations [deg]
theta = theta + tO;            % Add panel angle offset [deg]
theta = theta*(pi/180);        % Convert from degrees to radians [rad]

% Boundary points
X_b = cos(theta);              % Compute boundary point X-coordinate (radius of 1)
Y_b = sin(theta);              % Compute boundary point Y-coordinate (radius of 1)

% Number of panels
con_pt = length(X_b)-1;        % Number of panels (control points)

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points
edge = zeros(con_pt,1);        % Initialize con_pt
for i = 1:1:con_pt             % Loop over all panels
    edge(i) = (X_b(i+1)-X_b(i))*(Y_b(i+1)+Y_b(i)); % Compute edge value
end
sumEdge = sum(edge);           % Sum of all edge values
% If panels are CounterClockWise, flip them (don't if ClockWise)
if (sumEdge < 0)               % If panels are CCW
    X_b = flipud(X_b);           % Flip the X-data array
    Y_b = flipud(Y_b);           % Flip the Y-data array
end

%% PANEL METHOD GEOMETRY - REF [1]

% Initialize variables
X_cp   = zeros(con_pt,1);     % Initialize control point X-coordinate array
Y_cp   = zeros(con_pt,1);     % Initialize control point Y-coordinate array
L    = zeros(con_pt,1);       % Initialize panel length array
phiD = zeros(con_pt,1);       % Initialize panel orientation angle array [deg]

% Find geometric quantities of cylinder
for i = 1:1:con_pt                     % Loop over all panels
    X_cp(i)   = 0.5*(X_b(i)+X_b(i+1)); % X-value of control point
    Y_cp(i)   = 0.5*(Y_b(i)+Y_b(i+1)); % Y-value of control point
    dx      = X_b(i+1)-X_b(i);         % Change in X between boundary points
    dy      = Y_b(i+1)-Y_b(i);         % Change in Y between boundary points
    L(i)    = (dx^2 + dy^2)^0.5;       % Length of the panel
	phiD(i) = atan2d(dy,dx);           % Angle of the panel (positive X-axis to inside face)
    if (phiD(i) < 0)                   % Make all panel angles positive
        phiD(i) = phiD(i) + 360;
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD            =phiD+90;           % Angle from positive X-axis to outward normal vector [deg]
betaD             =deltaD-AoA;        % Angle between freestream vector and outward normal vector [deg]
betaD(betaD>360)=betaD(betaD>360)-360; % Make sure angles aren't greater than 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                 % Convert from [deg] to [rad]
beta = betaD.*(pi/180);                % Convert from [deg] to [rad]

%% COMPUTE SOURCE PANEL STRENGTHS - REF [5]

% Geometric integral (normal [I] and tangential [J])
% - Refs [2] and [3]
[I,J] = COMPUTE_SDM(X_cp,Y_cp,X_b,Y_b,phi,L); % Compute geometric integrals

% Populate A matrix
% - Simpler option: A = I + pi*eye(con_pt,con_pt);
A = zeros(con_pt,con_pt);     % Initialize the A matrix
for i = 1:1:con_pt            % Loop over all i panels
    for j = 1:1:con_pt        % Loop over all j panels
        if (i == j)           % If the panels are the same
            A(i,j) = pi;      % Set A equal to pi
        else                  % If panels are not the same
            A(i,j) = I(i,j);  % Set A equal to geometric integral
        end
    end
end

% Populate b array
% - Simpler option: b = -Uo*2*pi*cos(beta);
b = zeros(con_pt,1);         % Initialize the b array
for i = 1:1:con_pt           % Loop over all panels
    b(i) = -Uo*2*pi*cos(beta(i)); % Compute RHS array
end

% Compute source panel strengths (lambda) from system of equations
lambda  = A\b;               % Compute all source strength values

% Check the sum of the source strenghts
% - This should be very close to zero for a closed polygon
fprintf('Sum of L: %g\n',sum(lambda.*L));  

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities
% - Simpler method: Vt = Uo*sin(beta) + J*lambda/(2*pi);
%  Coefficient of pressure: Coeff_pres = 1-(Vt/Uo).^2;
Vt = zeros(con_pt,1);          % Initialize tangential velocity array
Coeff_pres= zeros(con_pt,1);   % Initialize pressure coefficient array
for i = 1:1:con_pt             % Loop over all i panels
    addVal  = 0;               % Reset the summation value to zero
    for j = 1:1:con_pt         % Loop over all j panels
        addVal = addVal + (lambda(j)/(2*pi))*(J(i,j)); % Sum all tangential source panel terms
    end
    
    Vt(i) = Uo*sin(beta(i)) + addVal;  % Compute tangential velocity by adding uniform flow term
    Coeff_pres(i) = 1-(Vt(i)/Uo)^2;    % Compute pressure coefficient
end

% Analytical angles and pressure coefficients
analyticTheta = linspace(0,2*pi,200)';             % Analytical theta angles [rad]
analyticCoeff_pres   = 1-4*sin(analyticTheta).^2;  % Analytical pressure coefficient []

%% COMPUTE LIFT AND DRAG

% Compute normal and axial force coefficients
CN = -Coeff_pres.*L.*sin(beta);                     % Normal force coefficient []
CA = -Coeff_pres.*L.*cos(beta);                     % Axial force coefficient []

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));       % Decompose axial and normal to lift coefficient []
CD = sum(CN.*sind(AoA)) + sum(CA.*cosd(AoA));       % Decompose axial and normal to drag coefficient []

% Display lift and drag coefficients in command window
fprintf('CL      : %g\n',CL);                       % Display lift coefficient (should be zero)
fprintf('CD      : %g\n',CD);                       % Display drag coefficient (should be zero)

%% COMPUTE STREAMLINES - REF [4]

if (flagPlot(4) == 1 || flagPlot(5) == 1)
    % Grid parameters
    nGridX = 300;                                   % X-grid for streamlines and contours
    nGridY = 300;                                   % Y-grid for streamlines and contours
    xVals  = [-1.5; 1.5];                           % X-grid extents [min, max]
    yVals  = [-1.5; 1.5];                           % Y-grid extents [min, max]

    % Streamline parameters
    stepsize = 0.01;                              % Step size for streamline propagation
    maxVert  = nGridX*nGridY*10;                   % Maximum vertices
    slPct    = 30;                                 % Percentage of streamlines of the grid
    Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';  % Create array of Y streamline starting points

    % Generate the grid points
    Xgrid   = linspace(xVals(1),xVals(2),nGridX)';  % X-values in evenly spaced grid
    Ygrid   = linspace(yVals(1),yVals(2),nGridY)';  % Y-values in evenly spaced grid
    [XX,YY] = meshgrid(Xgrid,Ygrid);                % Create meshgrid from X and Y grid arrays

    % Initialize velocities
    Vx = zeros(nGridX,nGridY);                      % Initialize X velocity matrix
    Vy = zeros(nGridX,nGridY);                      % Initialize Y velocity matrix

    % Solve for grid point X and Y velocities
    for m = 1:1:nGridX                              % Loop over X-grid points
        for n = 1:1:nGridY                          % Loop over Y-grid points
            XP = XX(m,n);                           % Isolate X point
            YP = YY(m,n);                           % Isolate Y point
            [Mx,My] = STREAMLINE_SPM(XP,YP,X_b,Y_b,phi,L); % Compute streamline Mx and My values (Ref [4])

            % Check if grid points are in object
            % - If they are, assign a velocity of zero
            [in,on] = inpolygon(XP,YP,X_b,Y_b);     % Find whether (XP,YP) is in polygon body
            if (in == 1 || on == 1)                 % If (XP, YP) is in the polygon body
                Vx(m,n) = 0;                        % X-velocity is zero
                Vy(m,n) = 0;                        % Y-velocity is zero
            else                                    % If (XP,YP) is not in the polygon body
                Vx(m,n) = Uo*cosd(AoA) + sum(lambda.*Mx./(2*pi));  % Compute X-velocity
                Vy(m,n) = Uo*sind(AoA) + sum(lambda.*My./(2*pi));  % Compute Y-velocity
            end
        end
    end
    % Compute grid point velocity magnitude and pressure coefficient
    Vxy  = sqrt(Vx.^2 + Vy.^2);         % Compute magnitude of velocity vector
    Coeff_presXY = 1-(Vxy./Uo).^2;      % Pressure coefficient []
end

%% PLOTTING

% FIGURE: Shape polygon with panel normal vectors
if (flagPlot(1) == 1)
    figure(1);                         % Create figure
    cla; hold on; grid off;            % Get ready for plotting
    set(gcf,'Color','White');          % Set color to white
    set(gca,'FontSize',12);            % Set font size
    fill(X_b,Y_b,'k');                 % Plot polygon
    for i = 1:1:con_pt                 % Loop over all panels
        X(1) = X_cp(i);                % Set X start of panel orientation vector
        X(2) = X_cp(i) + L(i)*cosd(betaD(i)+AoA);  % Set X end of panel orientation vector
        Y(1) = Y_cp(i);                            % Set Y start of panel orientation vector
        Y(2) = Y_cp(i) + L(i)*sind(betaD(i)+AoA);  % Set Y end of panel orientation vector
        plot(X,Y,'r-','LineWidth',3);              % Plot panel normal vector
    end
    xlabel('X Units');                              % Set X-label
    ylabel('Y Units');                              % Set Y-label
    axis equal;                                     % Set axes equal
    zoom reset;                                     % Reset zoom
end

% FIGURE: Geometry with the following indicated:
% - Boundary points, control points, first panel, second panel
if (flagPlot(2) == 1)
    figure(2);                                                   % Create figure
    cla; hold on; grid on;                                       % Get ready for plotting
    set(gcf,'Color','White');                                    % Set color to white
    set(gca,'FontSize',12);                                      % Set font size
    plot(X_b,Y_b,'k-','LineWidth',3);                             % Plot panels
    p1 = plot([X_b(1) X_b(2)],[Y_b(1) Y_b(2)],'g-','LineWidth',3); % Plot first panel
    p2 = plot([X_b(2) X_b(3)],[Y_b(2) Y_b(3)],'m-','LineWidth',3); % Plot second panel
    pB = plot(X_b,Y_b,'ko','MarkerFaceColor','k','MarkerSize',10); % Plot boundary points
    pC = plot(X_cp,Y_cp,'ko','MarkerFaceColor','r','MarkerSize',10); % Plot control points
    legend([pB,pC,p1,p2],...                                         % Show legend
           {'Boundary','Control','First Panel','Second Panel'});
    ylim([-1 1]);                                                   % Set Y-limits
    axis equal;                                                     % Set axes equal
    xlabel('X Units');                                              % Set X-label
    ylabel('Y Units');                                              % Set Y-label
    zoom reset;                                                     % Reset zoom
end

% FIGURE: Analytical and SPM pressure coefficient plot
if (flagPlot(3) == 1)
    figure(3);                                                      % Create figure
    cla; hold on; grid on;                                          % Get ready for plotting
    set(gcf,'Color','White');                                       % Set color to white
    set(gca,'FontSize',12);                                         % Set font size
    pA = plot(analyticTheta,analyticCoeff_pres,'k-','LineWidth',3); % Plot analytical pressure coefficient
    pC = plot(beta,Coeff_pres,'ks','MarkerFaceColor','r','MarkerSize',10); % Plot compute pressure coefficient
    xlabel('Angle [rad]');                                                 % Set X-label
    ylabel('Coeff_pres');                                                  % Set Y-label
    xlim([0 2*pi]);                                                        % Set X-limits
    ylim([-3.5 1.5]);                                                      % Set Y-limits
    legend([pA,pC],{'Analytical','SPM'},'Location','S');                    % Add legend
end

% FIGURE: Streamlines (and quiver if commented in)
if (flagPlot(4) == 1)
    figure(4);                                                              % Create figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
%     quiver(XX,YY,Vx,Vy,'r');                                              % Plot velocity vectors
    for i = 1:2:length(Ysl)                                                 % Loop over all streamlines
        sl = streamline(XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);    % streamline(X,Y,U,V,startx,starty)
        set(sl,'LineWidth',1);                                              % Change streamline width
    end
    fill(X_b,Y_b,'k');                                                      % Plot polygon
    xlim(xVals);                                                            % Set X-limits
    ylim(yVals);                                                            % Set Y-limits
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient contours
if (flagPlot(5) == 1)
    figure(5);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    contourf(XX,YY,Coeff_presXY,100,'EdgeColor','none');                    % Plot contour
    fill(X_b,Y_b,'k');                                                      % Plot polygon
    xlim(xVals);                                                            % Set X-limits
    ylim(yVals);                                                            % Set Y-limits
    zoom reset;                                                             % Reset zoom
end
