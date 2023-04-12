tic
hold on;

%----------system definition------------------
    T = 0.1;
    stepSize = 5
    vx = 5
    vy = 10

    A = [1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1]
    G = [T^2/2 0;T 0;0 T^2/2;0 T];
    C = [1 0 0 0;0 0 1 0]
    
    
    Qtilda = [0.3 0;0 0.1];
    zeroMeanProcessNoise = [0 0 0 0]
    zeroMeanMeasNoise = [0 0]
    
    Q = G*Qtilda*G'*1e6
    R = [0.1 0;0 0.1]*1e1;
    
    rng('default')  % For reproducibility
    wk = mvnrnd(zeroMeanProcessNoise,Q,stepSize)';
    vk = mvnrnd(zeroMeanMeasNoise,R,stepSize)';

%------------------system model------------------
    %acc = [wx;wy]
    %measNoise = [vx;vy]
    %xState = [x;vx;y;vy]
    %xStateNew = A*xState+G*acc
    %y = C*xState + measNoise

%----------True Position----------------------
    
    x = zeros(4,stepSize); %initialization
    x(2,:) = vx             % constant velocity assumption
    x(4,:) = vy
    y = zeros(2,stepSize);
    
    for i=1:stepSize
    
        x(:,i+1) = A*x(:,i) + wk(:,i)
        y(:,i) = C*x(:,i) + vk(:,i)
    end



%------------KF-----------------------
    estimation = [0;0]

    xzgz = x(:,1)
    pzgz = eye(4)
    yzgz = y(:,1)

     for i=1:stepSize
        clf;
        hold on;
        xlim([-5 stepSize*vx+5])
        ylim([-5 stepSize*vy+5])

        xogz = A*xzgz  % Predicton Update
        pogz = A*pzgz*A' + Q
        
        error_ellipse(C*pzgz*C' +R, C*xzgz,0.99)
        error_ellipse(C*pogz*C' +R, C*xogz,0.99)

        sogz = C*pogz*C' + R  % Measurement Update
        k1 = pogz*C'*inv(sogz)
        yhat1 = C*xogz
        estimation = [estimation yhat1]
        
        xogo = xogz + k1*(y(:,i)-yhat1)
        pogo = pogz - k1*sogz*k1'

        error_ellipse(C*pogo*C' + R, C*xogo,0.99)
       
        % for the next step new iterations
        xzgz = xogo
        pzgz = pogo
        pause(1)
        
        grid minor
        legend('Initial Random Data','Prediction Update','Measurement Update')
        title('Kalman Filter Implementation for Constant Velocity Model Step',i)

     end

toc

