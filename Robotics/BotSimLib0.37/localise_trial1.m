clf;        %clears figures
clc;        %clears console
clear;      %clears workspace
axis equal; %keeps the x and y scale the same

map=[0,0;60,0;60,45;45,45;45,59;106,59;106,105;0,105]; %default map
%maps{1} = [0,0;60,0;60,45;45,45;45,59;106,59;106,105;0,105]; %default map
%maps{2} = [0,0;60,0;60,50;100,50;70,0;110,0;150,80;30,80;30,40;0,80]; %long map
%maps{3} = [-30,0;-30,40;30,40;30,60;5,60;45,90;85,60;60,60;60,40;120,40;120,60;95,60;135,90;175,60;150,60;150,40;210,40;210,60;185,60;225,90;265,60;240,60;240,40;300,40;300,0]; %repeated features


% botSim = BotSim(map,[0.01,0.005,0]);  %sets up a botSim object a map, and debug mode on.
botSim = BotSim(map,[0,0,0]);  %sets up a botSim object a map, and debug mode on.
botSim.drawMap();
drawnow;
botSim.randomPose(10); %puts the robot in a random position at least 10cm away from a wall
target = botSim.getRndPtInMap(10);  %gets random target.

tic %starts timer

%your localisation function is called here.

%returnedBot = localise(botSim,map,target); %Where the magic happens

%we insert localisation in here to make it easier to edit
%% Localisation code

modifiedMap = map; %you need to do this modification yourself
botSim.setMap(modifiedMap);

targetp=BotSim(modifiedMap);
targetp.setBotPos(target);

isolated=1;
while isolated==1
    %generate some random particles inside the map
    num =300; % number of particles
    particles(num,1) = BotSim; %how to set up a vector of objects
    for i = 1:num
        particles(i) = BotSim(modifiedMap);  %each particle should use the same map as the botSim object
        particles(i).randomPose(0); %spawn the particles in random locations
    end

    position=[];
    heading=[];
    clusnum=floor(size(map,1)/2)-1; %assumed number of clusters
    if clusnum<2
        clusnum=2;
    end
    
    for i=1:num
        position=[position;particles(i).getBotPos()];
        heading=[heading;particles(i).getBotAng()];
    end
    
    [clusindex, roompos]=kmeans(position,clusnum); %split particles into clusters and get rooms

    adjacent=zeros(clusnum);
    n=100;
    R=5;
    for j=1:clusnum
        for i=1:clusnum
            if i==j
                score=0;
            else
                %Find the midpoint and determine if it can "see" its parents
                
                midpoint=(roompos(i,:)+roompos(j,:))/2; 
                newparticle=BotSim(modifiedMap);
                newparticle.setBotPos(midpoint);
                C=roompos(i,:)-midpoint;
                Cnorm=C/norm(C);
                theta=atan(Cnorm(2)/Cnorm(1));
                if Cnorm(1)<0
                    theta=theta+pi();
                end
                newparticle.setBotAng(theta);
                pscan=newparticle.ultraScan();
                if pscan(1)>=norm(roompos(i,:)-midpoint)
                    score1=1;
                else
                    score1=0;
                end
                %I can just use pscan(3) but I was thinking of expanding
                %code further
                
                D=roompos(j,:)-midpoint;
                Dnorm=D/norm(D);
                thetad=acos(Dnorm(1));
                if Dnorm(1)<0
                    thetad=thetad+pi();
                end
                newparticle.setBotAng(thetad);
                pscan=newparticle.ultraScan();
                if pscan(1)>=norm(roompos(j,:)-midpoint)
                    score2=1;
                else
                    score2=0;
                end
                score=score1*score2;
            end
        adjacent(i,j)=score;
        end
    end
    A=sum(adjacent, 2);
    if find(A==0) %a "room" is isolated. increase clusters
        clusnum=clusnum+1;
    else
        isolated=0;
    end
    
end
roompos
adjacent

roomparts(num,1) = BotSim;
for i=1:size(roompos,1)
    roomparts(i)=BotSim(modifiedMap);
    roomparts(i).setBotPos(roompos(i,:));
end


%HOW DO I CONNECT THE ROOMS
% midpoints=zeros(clusnum,clusnum,2);
% passable=zeros(clusnum);
% for j=1:clusnum
%     for i=1:clusnum
%         if i==j
%             midpoint=[0,0];
%         else
%             midpoint=(roompos(i,:)+roompos(j,:))/2;
%         end
%         
%         midpoints(i,j,:)=[midpoint];
%     end
% end
% 
% n=100;
% R = 5;
% score=zeros(size(midpoints,1),1);
% newparticle=BotSim(modifiedMap);
% for i=1:size(midpoints,1)
%     t = 2*pi * rand(n,1);
%     r = R*sqrt(rand(n,1));
%     points = midpoints(i,:)+[r.*cos(t), r.*sin(t)];
%     for j=1:size(points,1)        
%         newparticle.setBotPos(points(j,:));
%         if newparticle.insideMap()
%             score(i)=score(i)+1;
%         end
%     end
% end
% 
% score=score/n;

[targetdist, targetroom]=min(pdist2(roompos,target)); %find out which room target is in

maxNumOfIterations = 30;
n = 0;
converged =0; %The filter has not converged yet
particleScan=[];
mparticleScan=[];
rotation=[];
while(converged == 0 && n < maxNumOfIterations) %%particle filter loop
    n = n+1; %increment the current number of iterations
    botScan = botSim.ultraScan(); %get a scan from the real robot.
    [maxbS,maxbI]=max(botScan);   %identify index with max 

    %% Write code for updating your particles scans
    for i =1:num %for all the particles. 
        particleScan(:,i) = particles(i).ultraScan();
        [maxpS,maxpI]=max(particleScan(:,i));
        rotation(i)=maxbI-maxpI; %rotation matrix
        mparticleScan(:,i)=circshift(particleScan(:,i),rotation(i)); %align maxes
    end
    
    %% Write code for scoring your particles    
    
    particleScore=sqrt(sum((mparticleScan-botScan).^2,1));
       %particleScore=vecnorm(mparticleScan-botScan); %does not work on school
%    matlab
%    
%    particleScore = [];
%     for i =1:num %for all the particles. 
%         A=(mparticleScan(:,i)-botScan);
%         %A=mparticleScan(:,i)./botScan;
%         particleScore(i) = sqrt(dot(A,A)); %norm A
%     end
%     
    
    %% Write code for resampling your particles
    mparticleScore=particleScore; %duplicate particle score
    kill_list=zeros(num,1);    
    finite_list=isfinite(particleScore);
    index = find(finite_list==0); %finds the indices which are 0
    mparticleScore(index) = []; %remove infinites
    
    
    kill_list=particleScore>mean(mparticleScore); %kill particles further than mean
    
%     for i =1:num %for all the particles.
%         %higher score is bad    
%         if particleScore(i)>mean(mparticleScore)
%             kill_list(i)=1;
%         end
%     end
    respawn=sum(kill_list); %number to be respawned
    index = find(kill_list==1); %finds the indices which are 1
    particleScore(index) = [];
    particles(index)=[]; %Kill Particles!
    rotation(index)=[]; %maybe consider putting particles and rotation together but different datas
    
    
    farthest=max(particleScore);
    prob=1-particleScore/farthest;
    prob=prob/sum(prob);
    cumprob=cumsum(prob);
    
    for i=1:respawn %all the particles to respawn
        index=find(rand<=cumprob,1);
        position=particles(index).getBotPos();
        heading=particles(index).getBotAng()+rotation(index)*pi/3+0.20*randn;
        newparticle = BotSim(modifiedMap);
        newparticle.setBotPos(position);
        newparticle.setBotAng(heading);
        particles=[particles;newparticle];
    end    
    
    
    %% Write code to check for convergence   
    
    position=[];
    heading=[];
    for i=1:num
        position=[position;particles(i).getBotPos()];
        heading=[heading;particles(i).getBotAng()];
    end
    
    %clust=kmeans(position,2);
    botPosEs=mean(position);
    botAngEs=mod(mean(heading),2*pi);
    botPos=botSim.getBotPos();
    botAng=mod(botSim.getBotAng(),2*pi);
    
    disp(n)
    disp('where I think I am ')
    disp(botPosEs)
    disp(botAngEs)
    disp('where I really am')
    disp(botPos)
    disp(botAng)
    variance=norm(var(position))
    derror=norm(botPosEs-botPos)       
    herror=(botAngEs-botAng)/botAng        
    
    %% Write code to take a percentage of your particles and respawn in randomised locations (important for robustness)	
    five=floor(0.05*num);
    for i=1:five %five percent
        index=floor(num*rand)+1;
        particles(index)=[]; %kills it
        newparticle = BotSim(modifiedMap);
        newparticle.randomPose(min([min(botScan),10]));
        particles=[particles;newparticle];
    end 
    
    %% Write code to decide how to move next
    
    Z=target-botPosEs;
    targetdist=norm(Z);
    targetdir=atan(Z(2)/Z(1))-mean(heading);
    
    if mod(n,3)==1
        turn = (maxbI-1)*pi/3;
        move = min([maxbS-10,50]); 
    else
        turn = min([targetdir,1]);
        move = 2;
    end
    
    botSim.turn(turn); %turn the real robot.  
    botSim.move(move); %move the real robot. These movements are recorded for marking 
    for i =1:num %for all the particles. 
        particles(i).turn(turn); %turn the particle in the same way as the real robot
        particles(i).move(move); %move the particle in the same way as the real robot
    end
    
    
    %% Drawing
    %only draw if you are in debug mode or it will be slow during marking
    if botSim.debug()
        hold off; %the drawMap() function will clear the drawing when hold is off
        botSim.drawMap(); %drawMap() turns hold back on again, so you can draw the bots
        botSim.drawBot(30,'g'); %draw robot with line length 30 and green
        for i =1:num
            particles(i).drawBot(3); %draw particle with line length 3 and default color
        end
        targetp.drawBot(10,'r');
        for i =1:size(roompos,1)
            roomparts(i).drawBot(7, 'm'); %draw particle with line length 3 and default color
        end
        drawnow;
    end
end


returnedBot = botSim;



resultsTime = toc %stops timer

%calculated how far away your robot is from the target.
resultsDis =  distance(target, returnedBot.getBotPos())