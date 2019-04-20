function f=TSP_ACO()
x=[71,25,73,75,22,38,59,61,60,79,12,75,64,18,29,54,58,50,54,66,15,87,38,56,82,23,93,65,59,37,70,22,65,51,31,29,2,12,7,15,93,58,44,40,58,1,50,36,67,42,84,96,29,72,27,35,2,14,85,84,38,44,42,29,22,81,34,5,43,20,80,17,66,64,50,11,69,10,51,13,54,20,93,62,14,78,6,77,33,75,24,92,61,67,33,72,64,30,87,78];
y=[85,12,68,93,50,72,54,7,93,50,98,3,33,68,97,54,57,24,99,1,32,74,63,15,59,60,93,65,98,12,66,97,19,29,4,46,23,86,32,100,40,30,62,38,97,22,22,93,10,70,67,91,100,59,7,13,44,29,93,95,19,85,16,99,82,40,57,45,29,78,94,28,61,4,20,66,59,87,71,68,23,20,62,58,70,21,46,35,90,76,17,73,50,64,95,49,8,47,86,46];
C=[x;y];%positions of cities
n=length(C);
D=dist(C);

NC_max=50;%maximum iterations
m=200;%number of ants
alpha=1.4;%parameter that indicates the importance of pheromone
beta=2.2;%parameter that indicates the importance of heuristic factor
Rho=0.15;%coefficient that indicates pheromone evaporation 
Q=1e6;% coefficient that indicates increase of intensity of pheromone
%route_best
%length_best

%%main
Eta=1./D;
Tau=ones(n);
Route=zeros(m,n);
NC=1;
Length=zeros(m,1);
R_best=zeros(NC_max,n);
L_best=zeros(NC_max,1);
L_mean=zeros(NC_max,1);
set(gcf,'unit','normalized','position',[0.1,0.2,0.8,0.7]);
while NC<=NC_max
	%put m ants in n cities randomly
	pos=repmat(1:n,ceil(m/n));
	pos_ants=pos(1:m);
	pos_ants=pos_ants(:);
	Route(:,1)=pos_ants;% route
	
	%m ants choose next city according to the p(i,j)
	P=Tau.^alpha.*Eta.^beta;
	for i=1:m
		for j=2:n
			visited=Route(i,1:j-1);
			to_visit=1:n;
			to_visit(visited)=[];
			
			Psel=P(visited(end),to_visit);
			Pcum=cumsum(Psel);
			Pcum=Pcum/Pcum(end);
			index=find(Pcum>rand,1);
			Route(i,j)=to_visit(index);
		end
	end
	
	%refresh pheromone:Tau
	if NC>=2
		Route(1,:)=R_best(NC-1,:);
	end
	Delta_Tau=zeros(n);
	for i=1:m
		route=Route(i,:);
		for j=1:n-1
			p1=route(j);
			p2=route(j+1);
			d=D(p1,p2);
			Delta_Tau(p1,p2)=Delta_Tau(p1,p2)+Q/d;
			Length(i)=Length(i)+d;
		end
		p1=route(n);
		p2=route(1);
		d=D(p1,p2);
		Delta_Tau(p1,p2)=Delta_Tau(p1,p2)+Q/d;
		Length(i)=Length(i)+d;
	end
	Tau=(1-Rho).*Tau+Delta_Tau;
	
	%find best
	[L_best(NC),index]=min(Length);
	R_best(NC,:)=Route(index,:);
	L_mean(NC)=mean(Length);
	
	%plot
	subplot(1,2,1);
	route=R_best(NC,:);
	route(n+1)=route(1);
	scatter(C(1,:),C(2,:));
	hold on
	plot(C(1,route),C(2,route),'b');
	title('TSP问题优化结果 ')
	hold off
	
	subplot(1,2,2);
	plot(L_best,'b-')
	hold on
	plot(L_mean,'r-.')
	title('平均距离和最短距离') 
	legend('最短距离','平均距离');
	hold off
	pause(0.0001);
	%disp
	disp(['Iteration=' num2str(NC) ';MinDist = ' num2str(L_best(NC)) ';MeanDist = ' num2str(L_mean(NC))]);
	
	Route=zeros(m,n);
	Length=zeros(m,1);
	NC=NC+1;
end
R_best(end,:)
L_best(end)