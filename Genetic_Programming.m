clc;
clear all;

global modmat
global data

format long
tic
filedir='.\function1.csv'; 
data = textread(filedir, '', 'delimiter', ',');
% data(1,:)=data(1,:)+0.00000000000001;
totalrun=1000000;

file=load('.\generate30.mat');
population=file.evamatrix;

[poprow,popcol]=size(population);
ii=1;
while ii<=poprow
    matnow=population(ii,1:popcol-1).';
    modmat=matnow;
    esty=policy(1);
    inidist=distdif(esty);    
    population(ii,popcol)=inidist;        
    ii=ii+1;       
end

population=sortrows(population,popcol,'ascend');

endtempshort=zeros(poprow,totalrun+1);
endtempshort(:,1)=population(:,end); 

GAdistahead = population(1,end);

rr = 0;

summaryshort=zeros(totalrun+2,2);
summaryshort(1,1)=rr;
summaryshort(1,2)=GAdistahead;

while rr<=totalrun 
    pickrows=4;%round2even(randperm(10,1));
    matbegin=population(1:pickrows,:);     
    GAmatrix=crossover(matbegin);
    population((poprow-pickrows+1):end,:)=GAmatrix;
    population=sortrows(population,popcol,'ascend');  
    
    GAdistnow=population(1,end);    
    if GAdistnow<GAdistahead %shortest distance
        shortdist=GAdistnow;
    else
        shortdist=GAdistahead;      
    end               
    rr=rr+1;  
    endtempshort(:,rr+1)=population(:,end); 
    summaryshort(rr+1,1)=rr;
    summaryshort(rr+1,2)=shortdist;         
    GAdistahead = shortdist;       
end

modmat=population(1,1:popcol-1).';
esty=policy(1);
besty=esty;
test=showpolicy(1)    

shortx=summaryshort(:,1);
shorty=summaryshort(:,2);
figure(1)%plot shortest distance
plot(shortx,shorty,'b')

figure(2);    
hold on;
estimate=plot(data(:,1),besty,'b');    
question=plot(data(:,1),data(:,2),'r');
% ylim([-0.8 0.8])
grid on;
xlabel('X');
ylabel('Y');
title('Symbolic Regression');
% legend([estimate,question],{'Estimate','Question'},'Location','best');
hold off;  

goth=2; %run 4 times
savefilename=[filedir(1:end-4) 'GAselGo' num2str(goth) '.mat'];
save(savefilename)

toc

function GAmatrix=crossover(evamatrix)
global modmat
[matr,matc]=size(evamatrix);
evacrossmat=zeros(matr,matc);
ii=1;
while ii<=matr    
    gene1=evamatrix(ii,1:matc-1);
    modmat=gene1;
    treeg1=sort(shownodes(1));    
    gene2=evamatrix(ii+1,1:matc-1);  
    modmat=gene2;
    treeg2=sort(shownodes(1));      
    comtree=intersect(treeg1, treeg2); 
    rowcolcom=size(comtree);
    tt=randperm(rowcolcom(2),1);
    n=comtree(tt); 
    if n==1
        cronodes=1:63;
    else
        cronodes=n2n2n1(n);
    end        
    temp11=gene1(cronodes);
    temp12=gene2(cronodes);   
    gene1(cronodes)=temp12;
    newgene1=gene1;    
    gene2(cronodes)=temp11;
    newgene2=gene2;       
    matnow1=newgene1.';    
    guess=randperm(3,1);
    if guess==1|| guess==2     
        modmat=matnow1;
        treeahead=sort(shownodes(1));                 
        matnow1=mutat(matnow1,treeahead);  
    end            
    modmat=matnow1;
    esty=policy(1);
    newgene1dist=distdif(esty);         
    evacrossmat(ii,1:matc-1)= matnow1.';
    evacrossmat(ii,matc)= newgene1dist;         
    matnow2=newgene2.';    
    guess=randperm(3,1);
    if guess==1|| guess==2 
        modmat=matnow2;
        treeahead=sort(shownodes(1));
        matnow2=mutat(matnow2,treeahead);    
    end    
    modmat=matnow2;
    esty=policy(1);
    newgene2dist=distdif(esty);    
    evacrossmat(ii+1,1:matc-1)= matnow2.';
    evacrossmat(ii+1,matc)=newgene2dist;              
    ii=ii+2;   
end    
GAmatrix=evacrossmat;
end

function nodenum=shownodes(i)
global modmat
if modmat(i)==10000 %jia   
    strnum1=shownodes(2*i);
    strnum2=shownodes(2*i+1);   
    nodenum=[i,strnum1,strnum2];                  
elseif modmat(i)==20000 %jian      
    strnum1=shownodes(2*i);
    strnum2=shownodes(2*i+1);   
    nodenum=[i,strnum1,strnum2];                                   
elseif modmat(i)==30000 %cheng   
    strnum1=shownodes(2*i);
    strnum2=shownodes(2*i+1);   
    nodenum=[i,strnum1,strnum2];                    
elseif modmat(i)==40000 %chu 
    strnum1=shownodes(2*i);
    strnum2=shownodes(2*i+1);   
    nodenum=[i,strnum1,strnum2];        
elseif modmat(i)==50000%sin
    strnum1=shownodes(2*i); 
    nodenum=[i,strnum1];        
elseif modmat(i)==60000%cos  
    strnum1=shownodes(2*i); 
    nodenum=[i,strnum1];   
elseif modmat(i)==99999%cos       
    nodenum=i;     
else 
    nodenum=i;    
end
end

function answer=policy(i)
    
global modmat
global data
xvalues=data(:,1);

if modmat(i)==10000 %jia
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1+num2;                        
elseif modmat(i)==20000 %jian
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1-num2;                                         
elseif modmat(i)==30000 %cheng
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1.*num2;                           
elseif modmat(i)==40000 %chu
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1./num2;           
elseif modmat(i)==50000%sin
    num1=policy(2*i);   
    answer=sin(num1);        
elseif modmat(i)==60000%cos
    num1=policy(2*i);   
    answer=cos(num1); 
elseif modmat(i)==99999%cos  
    answer=xvalues;        
else
    answer=modmat(i);
end
end
 
function gongshi=showpolicy(i)
    
global modmat

if modmat(i)==10000 %jia   
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '+' strnum2 ')'];                  
elseif modmat(i)==20000 %jian      
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '-' strnum2 ')'];                                   
elseif modmat(i)==30000 %cheng   
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '*' strnum2 ')'];                      
elseif modmat(i)==40000 %chu 
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '/' strnum2 ')'];         
elseif modmat(i)==50000%sin
    strnum1=showpolicy(2*i); 
    gongshi=['(sin' strnum1 ')'];         
elseif modmat(i)==60000%cos  
    strnum1=showpolicy(2*i); 
    gongshi=['(cos' strnum1 ')'];    
elseif modmat(i)==99999%cos       
    gongshi='X';       
else
    answer=modmat(i);
    gongshi=['(' num2str(answer) ')'];    
end
end

function doutput=mutat(input63,tree)
    rowcol=size(tree);
    tt=randperm(rowcol(2),1);
    n=tree(tt);    
    if n==1 % ceng1=1;       
        doutput=new63mat();    
    else
        allsub=n2n2n1(n);
        for ii=allsub 
            if ii>=32
                doutput=putnu(input63,ii);
            else
                doutput=putfunu(input63,ii);                
            end
            input63=doutput;
        end
    end
end

function nn=n2n2n1(nodeset)
nodeahead=nodeset;
kk=1;%donot put 1st node
while kk<100       
    popnode=subnode(nodeahead);    
    if kk>=2^4 || popnode(1)>63
        break
    end        
    nodeset=[nodeset,popnode];
    kk=kk+1;       
    nodeahead=nodeset(kk);  
end
nn=nodeset;
end

function numhold=subnode(n)
    n1=2*n;
    n2=2*n+1;
    numhold=[n1,n2]; 
end

function routput63=new63mat()
    matsize=63;
    fuhao=[10000,20000,30000,40000,50000,60000];
    randmat=zeros(matsize,1);        
    for ii=1:31
        guess=randperm(2,1);
        if guess==1
            randmat(ii)=rand(1)*20-10;        
        elseif guess==2
            randmat(ii)=fuhao(randperm(6,1));                    
        end
    end
    for ii=32:matsize
        randmat(ii)=rand(1)*20-10;       
    end
    xhowmany=randperm(matsize,1);
    for jj=1:xhowmany
        xpick=randperm(matsize,1);
        randmat(xpick)=99999;
    end 
    routput63=randmat;    
end

function autput63=putfunu(fninput63,fnn)%Pick element to be any
    guess=randperm(2,1);
    if guess==1
        autput63=putfu(fninput63,fnn);        
    elseif guess==2
        autput63=putnu(fninput63,fnn);
    end
end    
        
function foutput63=putfu(finput63,fn)%Pick element to be sign
    % jia=10000;%string('+');
    % jian=20000;%string('-');
    % cheng=30000;%string('*');
    % chu=40000;%string('/');
    % san=50000;%string('sin');
    % kosan=60000;%string('cos');
    fuhao=[10000,20000,30000,40000,50000,60000];
    guess=randperm(6,1);
    if guess==1
        finput63(fn)=fuhao(1);        
    elseif guess==2
        finput63(fn)=fuhao(2);         
    elseif guess==3
        finput63(fn)=fuhao(3);         
    elseif guess==4
        finput63(fn)=fuhao(4);         
    elseif guess==5
        finput63(fn)=fuhao(5);         
    elseif guess==6
        finput63(fn)=fuhao(6);         
    end 
    foutput63=finput63;
end

function noutput63=putnu(ninput63,nn)%Pick element to be # or X
    guess=randperm(2,1);
    if guess==1
        ninput63(nn)=rand(1)*20-10;        
    elseif guess==2
        ninput63(nn)=99999; 
    end
    noutput63=ninput63;
end
        
function sumdist=distdif(inputy)
    global data
    yvalues=data(:,2);    
    [estyr,estyc]=size(inputy);
    if (estyr+estyc)>2
        sumdist=sum(abs(inputy-yvalues));
    else
        sumdist=Inf; 
    end
end 

function S = round2even(x) %cited from mathworks.com 
    if mod(x,2)<1 
        S = fix(x); 
    else 
        S =fix(x) + 1; 
    end
end
