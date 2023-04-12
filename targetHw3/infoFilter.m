tic
hold on;

%----------system definition------------------
    T = 0.1;
    stepSize = 5
    vx = 5
    vy = 10

    A = [1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1]
    G = [T^2/2 0;T 0;0 T^2/2;0 T];
    B=A;
    C = [1 0 0 0;0 0 1 0]
    
    
    Qtilda = [0.3 0;0 0.1];
    zeroMeanProcessNoise = [0 0 0 0]
    zeroMeanMeasNoise = [0 0]
    
    Q = G*Qtilda*G'*1e6;
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

    Zzgz =eye(4)% inv(pzgz)
    zzgz = inv(pzgz)*xzgz

     for i=1:stepSize
        clf;
        hold on;
%         xlim([-5 stepSize*vx+5])
%         ylim([-5 stepSize*vy+5])

        S = (inv(A)')*Zzgz*inv(A)
        J = S*B*inv(B'*S*B+inv(Q))
        K = Zzgz*C'*inv(R)


        zogz = (eye(4)- J*B')*inv(A)'*zzgz  % Predicton Update
        Zogz = (eye(4)- J*B')*S + Q
        
        pzgz = inv(Zzgz)
        xzgz = inv(Zzgz)*zzgz

        pogz = inv(Zogz)
        xogz = inv(Zogz)*zogz

        error_ellipse(C*pzgz*C' +R, C*xzgz,0.99)
        error_ellipse(C*pogz*C' +R, C*xogz,0.99)

% ---------------------Measurement Update------------------------
%         sogz = C*pogz*C' + R  
%         k1 = pogz*C'*inv(sogz)
         xogz = inv(Zogz)*zogz
         yhat1 = C*xogz
%         estimation = [estimation yhat1]
        
plot(y(1,i),y(2,i),'b*')
        zogo = zogz + C'*inv(R)*y(:,i)
        Zogo = Zogz + C'*inv(R)*C

        pogo = inv(Zogo)
        xogo = inv(Zogo)*zogo
        error_ellipse(C*pogo*C' + R, C*xogo,0.99)
       
        % for the next step new iterations
        zzgz = zogo
        Zzgz = Zogo
        pause(1)
        
        grid minor
        legend('Initial Random Data','Prediction Update','Actual Position','Measurement Update')
        title('Information Filter Implementation for Constant Velocity Model Step',i)

     end
toc


