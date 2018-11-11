function [Err]=EvaluateSL(stepLength,trueRec,gradient,model,source,dims)
recording = zeros(dims.nt,length(dims.recPos),'single');
Err=0;
for s = 1:dims.ds:length(dims.srcPos)
             u=zeros(dims.ny,dims.nx); %  time step n
             up1=zeros(dims.ny,dims.nx);% time step n+1
             um1=zeros(dims.ny,dims.nx);% time step n-1
            for t = 1:dims.nt
                 
                %  Solve wave equation using test model update
                [up1]=solveWaveEqn(u,up1,um1,t,dims,model+stepLength.*gradient,dims.srcPos(s),source);
                um1=u;
                u=up1;
                %  Record traces
                recording(t,:) = u(dims.recPos);
                
            end
            
            %% Calculate new error and check against old
            chi=recording-trueRec(:,:,s); %Residual
            phi=norm(chi);
            Err=sum(phi)+Err;
end           
            
            
        