clc;
clear all;

global matrix
global xvalues

format long
tic
filedir='.\function1.csv'; 
data = textread(filedir, '', 'delimiter', ',');
[rows,columns]=size(data);
xvalues=data(:,1);
yvalues=data(:,2);

% xvalues=[1,1;2,2];
totalrun=1000000;

jia=10000;%string('+');
jian=20000;%string('-');
cheng=30000;%string('*');
chu=40000;%string('/');
san=50000;%string('sin');
kosan=60000;%string('cos');

ceng1=1;
ceng2=2:3;
ceng3=4:7;
ceng4=8:15;
ceng5=16:31;
ceng6=32:63;

fuhao=[jia, jian, cheng, chu, san, kosan];
matrix=zeros(ceng6(end),1);

rr = 0;
bestyr=0;
goodcount=1;
closedist=0;
distahead=Inf;
summarylong=zeros(totalrun+1,2);

best100=zeros(ceng6(end),totalrun+1);

while rr<=totalrun   
    for ii=[ceng1 ceng2 ceng3 ceng4 ceng5]
        guess=randperm(2,1);
        if guess==1
            matrix(ii)=rand(1)*20-10;        
        elseif guess==2
            matrix(ii)=fuhao(randperm(6,1));                    
        end
    end
    for ii=ceng6
        matrix(ii)=rand(1)*20-10;       
    end

    [rowmat,colmat]=size(matrix);
    xhowmany=randperm(ceng6(end),1);
    for jj=1:xhowmany
        xpick=randperm(rowmat,1);
        matrix(xpick)=99999;
    end
    esty=policy(1);
    [estyr,estyc]=size(esty);
    if (estyr+estyc)>2
        distnow=sum(abs(esty-yvalues));
        if distnow<distahead
            targetset=matrix;
            closedist=distnow;
            besty=esty;            
            [bestyr,bestyc]=size(besty);
            best100(:,goodcount)=best100(:,goodcount)+matrix;
            goodcount=goodcount+1; 
        else
            closedist=distahead;                          
        end
        distahead = closedist;        
    end
    
    summarylong(rr+1,1)=summarylong(rr+1,1)+rr;
    summarylong(rr+1,2)=summarylong(rr+1,2)+closedist;    
    rr=rr+1;           
end

figure(1)
longx=summarylong(:,1);
longy=summarylong(:,2);
plot(longx,longy,'b')

if bestyr>2  
    figure(2);    
    hold on;
    estimate=plot(xvalues,besty,'b');    
    question=plot(xvalues,data(:,2),'r');
    % ylim([-0.8 0.8])
    grid on;
    xlabel('X');
    ylabel('Y');
    title('Symbolic Regression');
    legend([estimate,question],{'Estimate','Question'},'Location','best');
    hold off;   
end

matrix=targetset;
test=showpolicy(1)

goth=4; %run 4 times
savefilename=[filedir(1:end-4) 'randomGo' num2str(goth) '.mat'];
save(savefilename)

toc

function answer=policy(i)
    
global matrix
global xvalues

if matrix(i)==10000 %jia
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1+num2;                        
elseif matrix(i)==20000 %jian
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1-num2;                                         
elseif matrix(i)==30000 %cheng
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1.*num2;                           
elseif matrix(i)==40000 %chu
    num1=policy(2*i);
    num2=policy(2*i+1);
    answer=num1./num2;           
elseif matrix(i)==50000%sin
    num1=policy(2*i);   
    answer=sin(num1);        
elseif matrix(i)==60000%cos
    num1=policy(2*i);   
    answer=cos(num1); 
elseif matrix(i)==99999%cos  
    answer=xvalues;        
else
    answer=matrix(i);
end
end
 
function gongshi=showpolicy(i)
    
global matrix

if matrix(i)==10000 %jia   
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '+' strnum2 ')'];                  
elseif matrix(i)==20000 %jian      
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '-' strnum2 ')'];                                   
elseif matrix(i)==30000 %cheng   
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '*' strnum2 ')'];                      
elseif matrix(i)==40000 %chu 
    strnum1=showpolicy(2*i);
    strnum2=showpolicy(2*i+1);   
    gongshi=['(' strnum1 '/' strnum2 ')'];         
elseif matrix(i)==50000%sin
    strnum1=showpolicy(2*i); 
    gongshi=['(sin' strnum1 ')'];         
elseif matrix(i)==60000%cos  
    strnum1=showpolicy(2*i); 
    gongshi=['(cos' strnum1 ')'];    
elseif matrix(i)==99999%cos       
    gongshi='X';       
else
    answer=matrix(i);
    gongshi=['(' num2str(answer) ')'];    
end
end

% % eval('cos((3,141)()')
