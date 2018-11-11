%% Setting up dimensions
dims.dy =     10; % [m]
dims.dx =     10; % [m]
dims.dt = 1.0e-3; % [s]

dims.ny = 201; % Cells in y-direction
dims.nx = 301; % Cells in x-direction
dims.nt = 801; % Amount of time steps 

%% Model dimensions
dims.modely = 100:150;
dims.modelx = 100:200;
dims.my = length(dims.modely);
dims.mx = length(dims.modelx);

%% Source locations
sx = min(dims.modelx):max(dims.modelx);
sy = min(dims.modely)*ones(1,length(sx));
dims.srcPos = sy + dims.ny*sx;

%% Receiver locations
rx = min(dims.modelx):max(dims.modelx);
ry = min(dims.modely)*ones(1,length(rx));
dims.recPos = ry+dims.ny*rx;

%% Creating background model
bg = zeros(dims.ny,dims.nx,'single');
bg(:) = 2.0e3;         % [m/s] - Background
bg(115:end,:) = 2.3e3; % [m/s] - Layer

%% Begin iteration
model = bg;     % Starting model
dims.ds = 10;   % Grid point distance between sources
maxIter = 10;   % Maximum number of iterations per frequency
freqs = [12];  % Frequencies to use in inversion

errVec = zeros(1,maxIter*length(freqs));
alphas=zeros(1,maxIter*length(freqs));

it = 1; tic;
for f = freqs
    %% Generating ricker source signature wavelet 
    source = rickerWave(f,dims);
    %plot(1:801,source);
    %% Load true recording
    load (['trueRec_',num2str(f),'Hz.mat']);   
    for i = 1:maxIter
        %% Calculate gradient ## IMPLEMENT ##
        [gradient,err] = calculateGradient(trueRec,source,model,dims) ;
        
        figure()
        subplot(121)
        imagesc(gradient(dims.modely,dims.modelx))
        title('Gradient')
        xlabel('Distance')
        ylabel('Depth')
        colorbar
        drawnow()

            
        %% Taper gradient
        gradient = taperGradient(gradient);
        subplot(122)
        imagesc(gradient(dims.modely,dims.modelx))
        title('Gradient')
        xlabel('Distance')
        ylabel('Depth')
        colorbar
        drawnow()
        

        %% Calculate step length ## IMPLEMENT ##
        [stepLength,err] = calculateStepLength2(model,trueRec,source,dims,gradient,err);
        
        alphas(it)=stepLength;
              
        
        %% Update model
        model = model + stepLength*gradient;

        errVec(it) = err; 
        
        
        figure(3)
        imagesc(model(dims.modely,dims.modelx))
        title('Model Updated')
        xlabel('Distance')
        ylabel('Depth')
        colorbar
        drawnow()
        
        figure(4)
        subplot(2,1,1)
        semilogy(errVec,'-')
        xlabel('Iteration')
        ylabel('Objective function')
        drawnow()
        subplot(2,1,2)
        semilogy(alphas,'-')
        xlabel('Iterations')
        ylabel('Step Length') 
        drawnow();
               
        it = it + 1;
                      
        toc
       
        
        
    end
    
end


        
