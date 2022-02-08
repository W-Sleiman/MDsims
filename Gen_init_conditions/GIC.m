%clear all
%close all
%clc

numRun = 100;
runID = 1;

while ( runID <= numRun)
    
    N = 1000;
    Coordinates = zeros(N,3);
    BigMatrix = false(N,3*N/2,3*N/2);
    BigMatrix(N,3*N/4,3*N/4) = 1;
    TempCoord = [N,3*N/4,3*N/4];
    
    i = 1;
    failSafe = 0;
    tol = 1e4;
    checker = 1;

    while i <= N

        XYZ = ceil(3*rand);
        step = round(rand);
        if (step==0)
            step=-1;
        end
        TempCoord(XYZ) = TempCoord(XYZ)-step;
        
        if (TempCoord(1)<N)
            if(BigMatrix(TempCoord(1),TempCoord(2),TempCoord(3))==0)
                Coordinates(i,:) = TempCoord;
                BigMatrix(TempCoord(1),TempCoord(2),TempCoord(3)) = 1;
                i = i+1;
            else
                TempCoord(XYZ) = TempCoord(XYZ)+step;
            end
        else
            TempCoord(XYZ) = TempCoord(XYZ)+step;
        end

        failSafe = failSafe + 1;
        
        if (failSafe==tol)
            i = N+1;
            checker = 0;
        end
        
    end

    runID

    if (checker == 1)

        BeforeFinalCoord = Coordinates(Coordinates(:,1)~=0,:);
        FinalCoord = zeros(N,3);
        FinalCoord(:,1) = (BeforeFinalCoord(1:N,2)-750);
        FinalCoord(:,2) = BeforeFinalCoord(1:N,3)-750;
        FinalCoord(:,3) = -1*(BeforeFinalCoord(1:N,1)-1000);
        FinalCoord(:,4) = 2*rand(N,1)-1;
        FinalCoord(:,5) = 2*rand(N,1)-1;
        FinalCoord(:,6) = 2*rand(N,1)-1;
        
        for j = 1:N-1
            FinalCoord(N-j+1,:)=FinalCoord(N-j,:);
        end

        FinalCoord(1,:) = 0;

        csvwrite(strcat('InitialConditions/P_start_',num2str(runID)),FinalCoord);

        A = FinalCoord;
%         figure(999);plot3(A(:,1),A(:,2),A(:,3),'mo',A(:,1),A(:,2),A(:,3),'b')

        runID = runID + 1;

    end

end