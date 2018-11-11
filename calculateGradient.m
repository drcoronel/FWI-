function [gradient, err] = calculateGradient(trueRec,source,model,dims) 
    recording = zeros(dims.nt,length(dims.recPos),'single');
    gradient  = zeros(dims.ny,dims.nx,'single');
    forwardField = zeros(dims.my,dims.mx,dims.nt,'single'); 
    adjointField = zeros(dims.my,dims.mx,dims.nt,'single');
    rec=zeros(dims.nt,length(dims.recPos),length(dims.srcPos));
    err = 0;
    for s = 1:dims.ds:length(dims.srcPos)
        %% Run forward simulation on background model
        u=zeros(dims.ny,dims.nx); %  time step n
        uap1=zeros(dims.ny,dims.nx);% time step n+1
        um1=zeros(dims.ny,dims.nx);% time step n-1
        for t = 1:dims.nt
            % Solve wave equation'
            [uap1]=solveWaveEqn(u,uap1,um1,t,dims,model,dims.srcPos(s),source);
            um1=u;
            u=uap1;
            % Record traces
            recording(t,:) = u(dims.recPos);
            % Save fprward field for use in correlation
            forwardField(:,:,t) = u(dims.modely,dims.modelx); 
        end

        %% Calculate difference and error
        rec(:,:,s)=recording;% Save recordings for each source.
        chi=recording-trueRec(:,:,s); %Residual            
        phi=norm(chi);  %objective function
        err=sum(phi)+err;
        adchi=flipud(chi); %Inverse residual used for the adjoint simulation                      
       
        
        %% Run adjoint simulation
        ua=zeros(dims.ny,dims.nx); %  time step n for adjoint field 
        uap1=zeros(dims.ny,dims.nx);% time step n+1
        uam1=zeros(dims.ny,dims.nx);% time step n-1
        for t = 1:dims.nt
        
            % Solve wave equation using the difference (chi) as sources
            [uap1]=solveWaveEqn(ua,uap1,uam1,t,dims,model,dims.recPos,adchi);
            uam1=ua;
            ua=uap1;
            % Save adjoint field for use in correlation
            adjointField(:,:,dims.nt-t+1) = ua(dims.modely,dims.modelx);
        end
        %% Correlate
        for t = 2:dims.nt-1
            % Calculate the time derivative of the displacement to
            % gradient.
            fF2=forwardField(:,:,t+1);fF1=forwardField(:,:,t-1);
            dtfwd=(fF2-fF1)./(2*dims.dt); 
            
            aF2=adjointField(:,:,t+1);aF1=adjointField(:,:,t-1);
            dtadj=(aF2-aF1)./(2*dims.dt);
            
            gradient(dims.modely,dims.modelx)=gradient(dims.modely,dims.modelx)+(dtfwd.*dtadj);        
            
            
        end
    end
   
end

