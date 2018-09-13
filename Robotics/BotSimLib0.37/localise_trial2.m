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

scannum = 6;
botSim.setScanConfig(botSim.generateScanConfig(scannum)); 
cycle=0;

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
                %I can just use pscan(3) but changes with differing scannum
                
                D=roompos(j,:)-midpoint;
                Dnorm=D/norm(D);
                theta=atan(Dnorm(2)/Dnorm(1));
                if Dnorm(1)<0
                    theta=theta+pi();
                end
                newparticle.setBotAng(theta);
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
        num=num+50;
    else
        isolated=0;
    end
    [targetdist, targetroom]=min(pdist2(roompos,target)); %find out which room target is in
    E=target-roompos(targetroom);
    Enorm=E/norm(E);
    theta=atan(Enorm(2)/Enorm(1));
    if Enorm(1)<0
        theta=theta+pi();
    end
    targetp.setBotAng(theta);
    pscan=targetp.ultraScan();
    if pscan(1)<=norm(target-roompos(targetroom))
        %target is in the wrong room!
        clusnum=clusnum+1;
        num=num+50;
        isolated=1; %redo reseed
    end
end
roompos
adjacent

% just some dummy particles so I can see rooms
roomparts(num,1) = BotSim;
for i=1:size(roompos,1)
    roomparts(i)=BotSim(modifiedMap);
    roomparts(i).setBotPos(roompos(i,:)); 
end

% do wavefront/djvorak here
visited=[];
unvisited=1:size(roompos,1);
current=targetroom;
distancematrix=squareform(pdist(roompos));
djmatrix=zeros(size(roompos));%not really but they are same size
djmatrix(:,1)=Inf;
djmatrix(current,1)=0;

while ~isempty(unvisited)
    %djorak
    [dist, index]=min(djmatrix(unvisited,1)); %find minimum dist in unvisited
    current=unvisited(index);
    visitlist=find(adjacent(current,:)); %gives indexes adjacent to current
    [value, index]=intersect(visitlist,visited);
    visitlist(index)=[]; %kill off nodes already visited
    totaldistance=dist+distancematrix(current,visitlist);
    
    for index=find(totaldistance'<djmatrix(visitlist,1));
        djmatrix(visitlist(index),1)=totaldistance(index);
        djmatrix(visitlist(index),2)=current;
    end
    
    visited=[visited;current]
    unvisited(unvisited==current)=[]   
end
    
    
    
    


%find shortest distance to targetroom


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
        particles(i).setScanConfig(particles(i).generateScanConfig(scannum)); 
        particleScan(:,i) = particles(i).ultraScan();
        [maxpS,maxpI]=max(particleScan(:,i));
        rotation(i)=maxbI-maxpI; %rotation matrix
        mparticleScan(:,i)=circshift(particleScan(:,i),rotation(i)); %align maxes
    end
    
    %% Write code for scoring your particles    
    
    particleScore=sqrt(sum((mparticleScan-botScan).^2,1));
    particleScore=particleScore.^-1;
%     %norm partscan-botscan
%     %is this weight?
%     
%     mparticleScore=particleScore; %duplicate particle score
%     
%     finite_list=isfinite(particleScore);
%     index = find(finite_list==0); %finds the indices which are 0
%     particleScore(index) = []; %remove infinites
%     particles(index)=[];%kill particles with inf. Outsiders!
%     rotation(index)=[]; %particles, its score, its rotation still the same
%     %prob=exp(-(particleScore).^2/(2*var(particleScore)^2));
%     %prob=0.5*exp(-(particleScore)/2);
%     
%     
%     farthest=max(particleScore);
%     prob=1-particleScore/farthest;
%     prob=prob.^2; %This is just to make convergence faster
%     prob=prob/sum(prob); %normalized this might be weight now
%     cumprob=cumsum(prob);
%     
%     
%     corr=zeros(num,1);
%     particleScore=zeros(num,1);
%     for i=1:num
%         C=cov(botScan,mparticleScan(:,i));
%         corr(i)=C(1,2)/sqrt(C(1,1)*C(2,2));
%     end
%     particleScore=10*exp(-(corr-1).^2);
%     finite_list=isfinite(particleScore);
%     index = find(finite_list==0); %finds the indices which are 0
%     particleScore(index) = []; %remove infinites
%     particles(index)=[];%kill particles with inf. Outsiders!
%     rotation(index)=[]; %particles, its score, its rotation still the same
%     %prob=exp(-(particleScore).^2/(2*var(particleScore)^2));
%     %prob=0.5*exp(-(particleScore)/2);
%     
    prob=particleScore/sum(particleScore);
    index = find(prob==0); %finds the indices which are 0
    particleScore(index) = []; %remove infinites
    particles(index)=[];%kill particles with inf. Outsiders!
    rotation(index)=[]; %particles, its score, its rotation still the same
    prob(index)=[];
    cumprob=cumsum(prob);

    
    
    %% Write code for resampling your particles
    weight=zeros(num,1)+1;
    newset(num,1)=BotSim;
    U=1./num;
    B=rand*2*max(prob);
    if sum(cumprob)>0
        for i=1:num %all the particles to respawn
            index=find(B<=cumprob,1);
            B=B+U;
            if B>1
                B=B-1;
            end
            %index=find(rand<=cumprob,1);
            position=particles(index).getBotPos()+0.1*rand(1,2);
            heading=particles(index).getBotAng()+rotation(index)*2*pi/scannum+0.05*randn;
            newset(i) = BotSim(modifiedMap);
            newset(i).setBotPos(position);
            newset(i).setBotAng(heading);
            weight(i)=prob(index);
        end
    else
        for i=1:num %all the particles to respawn
            newset(i) = BotSim(modifiedMap);
            newset(i).randomPose(5);
            weight(i)=1;
        end
    end
    particles=newset;
    weight=weight/sum(weight); %returns prob total 1
    
    %% Write code to check for convergence   
    
    position=[];
    heading=[];
    for i=1:num
        position=[position;particles(i).getBotPos()];
        heading=[heading;particles(i).getBotAng()];
    end
    
    %clust=kmeans(position,2);
    botPosEs=sum(position.*weight);
    botAngEs=mod(sum(heading.*weight),2*pi);
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
    distrib=norm(botPosEs-mean(position)) %this can work as convergence!
    
    
    
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
    
    confident=0;
    if distrib<1
        confident =1; % the bot knows where it is! I hope
    end
    
    if confident==0
        %robot is lost
        %move randomly to localize
        if mod(cycle,3)==0
            turn = (maxbI-1)*2*pi/scannum;
            move = min([maxbS-10,50]); 
        else
            turn=0.5;
            move = 5;
        end
        cycle=cycle+1;
    else
        %robot kinda knows where it is
        %Lets try djestra or whatever
        %A* from there
        %nodes are roompos
        cycle=0;
        %find out robo room
        
        [robodist,roboroom]=pdist2(roompos,botPosEs,'euclidean','Smallest',1);
        current=roboroom;
        
        
        %head to nearest node
        if robodist>3
            travelto=roompos(roboroom,:);
        else
            if roboroom==targetroom
                travelto=target;
                converged=1;
            else
                travelto=roompos(djmatrix(current,2),:);
            end
        end
        
        %pairwise distance between nodes
        %squareform(pdist(roompos))./adjacent
        
        %euclid dist of rooms to target room
        %sqrt(sum((roompos-roompos(targetroom)).^2,2))
        %adjacentcy matrix is adjacent
        %find best path from roboroom to targetroom
        %once in targetroom go to target
        
        Z=travelto-botPosEs;
        Zdist=norm(Z);
        Zdir=atan(Z(2)/Z(1))-botAngEs;
        if Z(1)<0
            Zdir=Zdir+pi;
        end
        
        move=Zdist
        turn=Zdir
        
    end
    

    
    

    
    
    %collision detection
    collide=botScan<10;
    collision=1;
    while collision==1
        dir=turn/(2*pi/scannum);
        dir=round(dir+1);
        if dir<=0
            dir=dir+6
        end
        if collide(dir)
            turn=rand*2*pi-0.001;
            move=5;
        else
            collision=0;
        end
            
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
        
        %mostly because I do not know how to draw an x
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