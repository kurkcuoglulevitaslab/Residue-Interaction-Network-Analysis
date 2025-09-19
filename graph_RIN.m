%Code by Dr. Ozge Kurkcuoglu
%Istanbul Technical University, 2025
%The code requires at least MATLAB R2019a

clear all
clc


A=load('connectivity.out'); % create a graph based on contact topology
C=pdbread('coordmin1.pdb');
res=C.Model.Atom(end).AtomSerNo; %res=number of residues

for i=1:res
x(i,1)=C.Model.Atom(i).X;
y(i,1)=C.Model.Atom(i).Y;
z(i,1)=C.Model.Atom(i).Z;
end

s=A(:,1); %index
t=A(:,2); %local interaction strength aij
icount=0;

weights=1./A(:,3); %edge weight 

names=sprintfc('%g',1:res);

G=graph(s,t,weights,names); %construct the network
G.Nodes.XCoord=x;
G.Nodes.YCoord=y;
G.Nodes.ZCoord=z;

%show the network
p=plot(G,'XData',x,'YData',y,'ZData',z,'MarkerSize',5)

%calculate betweenness centrality and show the values on the plot
wbc=centrality(G,'betweenness','Cost',G.Edges.Weight);
n=numnodes(G);
p.NodeCData=2*wbc./((n-2)*(n-1));

colormap jet
colorbar
title('RIN colored with respect to CB values')
%---------------------------

resno=1:res;
yy=p.NodeCData;

figure(2)
plot(resno,yy)
xlabel('residue index')
ylabel('betweenness centrality')
title('CB values of the residues')

figure(3)
hist(yy,10)
xlabel('betweenness centrality')
ylabel('frequency')
title('CB distribution')

TOP5=quantile(yy,0.95) %threshold value to select hub residues

fileID = fopen('betweenness.out','w');
fileID2= fopen('top_5_percent.out','w');

fprintf(fileID,'i resn ch  CB     \n');
fprintf(fileID2,'i resn ch  CB     \n');

for i=1:res
    fprintf(fileID,'%d %s  %s  %9.7f\n',C.Model.Atom(i).resSeq,C.Model.Atom(i).resName,C.Model.Atom(i).chainID,yy(i));
    
    if yy(i) >= TOP5
        fprintf(fileID2,'%d  %s  %s  %9.7f\n',C.Model.Atom(i).resSeq,C.Model.Atom(i).resName,C.Model.Atom(i).chainID,yy(i));
    end
end

fclose(fileID);
fclose(fileID2);

%Spectral analysis to partition the network
L=laplacian(G);
 [V,D]=eigs(L,11,'sm');
 w=V(:,1);
 mode1= V(:,2);
 mode2= V(:,3);
 mode3= V(:,4);
 mode4= V(:,5);
 mode5= V(:,6);
 mode6= V(:,7);
 mode7= V(:,8);
 mode8= V(:,9);
 mode9= V(:,10);
 mode10= V(:,11);

 mod = zeros(res,10);

for i=1:10
    mode_n = eval(sprintf('mode%d',i));
    mod(:,i)=mode_n(:,1);
end

nmod = zeros(res,10);
% normalization for clarity
for i=1:10
    for j=1:res
        nmod(j,i)= mod(j,i)/sqrt(mod(j,i)*mod(j,i));
    end
end

% write values into an out file
writematrix(mod,'10modes.txt','Delimiter','tab')

% write normalized values into an out file
writematrix(nmod,'10modes_normalized.txt','Delimiter','tab')

%show the partitioning of the network at the smallest non-zero eigen:
%Fiedler vector, i.e. mode1
figure(4)
p=plot(G,'XData',x,'YData',y,'ZData',z,'MarkerSize',5)
highlight(p,find(mode1>=0),'NodeColor','r')
highlight(p,find(mode1<0),'NodeColor','k')
title('Spectral analysis, Fiedler vector')





