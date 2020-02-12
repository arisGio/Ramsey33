%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% START OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

rng("shuffle")            % randomness by instance
n=7;                      % # of vertices [ min = 2 ]
trials=1000000;           % # of trials #50 is good
beta = 2;                 % temperature in reduced units
B=ones(n)-eye(n);         % bonds matrix
first = three_energy(n,B);
energy = zeros(1,trials); % preallocation
for i=1:trials
    energy(i) = three_energy(n,B); 
    %%%%%%%%%%%%%%%%% randomd update (via indices) B %%%%%%%%%%%%%
    row = randi(n-1); % last row is not permitted
    col = randi([min(row+1,n),max(row+1,n)]); % looks obv. but it's not
    B(row,col) = B(row,col)*(-1);
    B(col,row) = B(row,col);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energy(i+1) = three_energy(n,B); 
    delta_energy = energy(i+1) - energy(i); % energy difference
    if rand < acceptance(beta, delta_energy)
        continue
    else
        B(row, col) = B(row, col)*(-1); % updating bonds matrix
        B(col,row) = B(row,col);
        energy(i+1) = energy(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
semilogx(1:length(energy),energy,'*r')
figure(2)
%%%%%%%%%%%%%%%%%%%%%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input graph weights: Bonds matrix will give us the weights
G=graph(B,'upper');
plot(G,'NodeColor','k','EdgeCData',G.Edges.Weight)
c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
colormap cool
%%%%%%%%%%%%%%%%%%%%%% END GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%

toc

%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = three_energy(n,B) % energy of all possible triangles
    e = 0;
    for i=1:n
        for j=(i+1):n
            for k=(j+1):n
                e = e + (B(i,j)+ B(j,k)+B(k,i))^2;
            end
        end
    end
end
             
function a = acceptance(beta, delta_energy) % acceptance rule
    a = min(1,exp(-beta*delta_energy));
end
%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
 
 
 
