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
    
    Q = G*Qtilda*G'
    R = [0.1 0;0 0.1];
    
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
%         xlim([-5-10 stepSize*T*vx+5+100])
%         ylim([-5-10 stepSize*T*vy+5+100])
        xlim([-5 11])
        ylim([-6 13])

        xogz = A*xzgz  % Predicton Update
        pogz = A*pzgz*A' + Q
        
        error_ellipse(C*pzgz*C' +R, C*xzgz,0.99)
        error_ellipse(C*pogz*C' +R, C*xogz,0.99)

        sogz = C*pogz*C' + R  % Measurement Update
        k1 = pogz*C'*inv(sogz)
        yhat1 = C*xogz
        estimation = [estimation yhat1]

        falseAlarmNum = 10
        sz = [1 falseAlarmNum];
        gateX = unifrnd(-5, stepSize*T*vx+5,sz)
        gateY = unifrnd(-5, stepSize*T*vy+5,sz)

       c = 0
       cY = []
      
          gate = [gateX(1,:); gateY(1,:)]
          minNorm = 1000
%         for k = 1:falseAlarmNum
%             if norm(gate(:,k)) < minNorm
%                 minNorm = norm(gate(:,k))
%                 c = k
%             end
%         end
        
%         for k = 1:stepSize
%             if (norm([gateX(1,k); gateY(1,k)]) == minNorm)
%                 c = k
%             end
% %             if (abs(gateY(1,k)) == min(abs(gateY(1,:))))
% %                 cY = k
% %             end
%         end
        if (i<=stepSize-1)
            plot(y(1,i+1),y(2,i+1),'b*')
        end
%         y(1,i) = gateX(1,c) + y(1,i)
%         y(2,i) = gateY(1,c) + y(2,i)
%         plot(y(1,i),y(2,i),'g*')

        gammaG= chi2inv(0.9,2)
gateMinimizer =[]
        for k = 1:falseAlarmNum
            if ((gate(:,k)-yhat1)'*inv(sogz)*(gate(:,k)-yhat1) < gammaG)
                    if  (norm(y(:,i) - gate(:,k))< minNorm)
                        minNorm = norm(y(:,i) - gate(:,k))
                        y(:,i) = gate(:,k)
%                         cY = [cY k]   
                    end     
            end
%             if ((y(:,i)+gate(:,k)-yhat1)'*inv(sogz)*(y(:,i)+gate(:,k)-yhat1) < gammaG)
%                 
% %                 if (norm([gateX(1,k); gateY(1,k)]) < minNorm)
%                     y(:,i) = y(:,i)+gate(:,k)
%                     cY = [cY k]   
%             % end     
%             end
        end
%         [M,index] = min(minNorm)
%        circle([y(1,i)- y(2,i)],minNorm,1000,'--m')

%         if length(cY)>1
%             for m = 1:length(cY)
%                 gateMinimizer = [gateMinimizer norm(gate(:,cY(m)))]
%                 [MinGate,I] = min(gateMinimizer)
%                 y(:,i) = y(:,i)+gate(:,cY(I))
%             end
%         end
        plot(y(1,i),y(2,i),'g*')

        xogo = xogz + k1*(y(:,i)-yhat1)
        pogo = pogz - k1*sogz*k1'

        error_ellipse(C*pogo*C' + R, C*xogo,0.99)
        %estimation = [estimation C*xogo]
        plot(estimation(1,:),estimation(2,:))

        for j = 1:falseAlarmNum
            plot(gateX(1,j)+y(1,i),gateY(1,j)+y(2,i),'r*')
        end

        % for the next step new iterations
        xzgz = xogo
        pzgz = pogo
        pause(1)
        
        grid minor
        if (i<=stepSize-1)
            legend('Initial Random Data','Prediction Update','Actual Measurement','Gating Result','Measurement Update')%,'Uniformly Distributed False Alarms')
        else
            legend('Initial Random Data','Prediction Update','Actual Measurement','Measurement Update')%,'Uniformly Distributed False Alarms')
        end
            title('Kalman Filter Implementation with Gating (Number of False Alarms =10) for Constant Velocity Model Step',i)

     end



