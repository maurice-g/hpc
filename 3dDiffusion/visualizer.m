function visualizer()
%%%%%%%%%% define some general parameters %%%%%%%%%%
file_basename = 'densities-';

dx=.25;

m=100;
n=100;
p=100;

topology = [2;2;2];

ml = m/topology(1);
nl = n/topology(2);
kl = p/topology(3);

densities = zeros(m,n,p);

%%%%%%%%%% loop through all files (one per core) and read content %%%%%%%%%
size_string = ['_' num2str(ml) '-' num2str(nl) '-' num2str(kl) '.txt'];
tmp = zeros(ml,nl,kl);
for i=0:(topology(1)-1)
    for j=0:(topology(2)-1)
        for k=0:(topology(3)-1)
            filename = [file_basename num2str(i) '-' num2str(j) '-' num2str(k) size_string];
            disp(['Reading file: ' filename]);
            data = load(filename);
            for depth=1:kl
                tmp(:,:,depth) = data((depth-1)*nl+1:depth*nl,:)';
            end
            densities(i*ml+1:(i+1)*ml,j*nl+1:(j+1)*nl,k*kl+1:(k+1)*kl) = tmp;
        end
    end
end
        
%%%%%%%%%% plot the resulting 3D-array %%%%%%%%%%
[x,y,z] = meshgrid(dx*(0:m-1),dx*(0:n-1),dx*(0:p-1));
%[x,y,z] = meshgrid(dx*(0:m-1));


%scatter3(x(:),y(:),z(:),40,densities(:), 'filled');
slice(x,y,z,densities,[.5;1.25],[.5],[1])


colorbar

uicontrol('Style', 'slider',...
        'Min',0,'Max',dx*(p-1),'Value',.25,...
        'Position', [400 20 120 20],...
        'Callback', {@surfzlim,x,y,z,densities});   
    
    
end

function surfzlim(hObj,event,x,y,z,densities)
    % Called to set zlim of surface in figure axes
    % when user moves the slider control
    val = get(hObj,'Value');
    slice(x,y,z,densities,[10],[10],val);
end






