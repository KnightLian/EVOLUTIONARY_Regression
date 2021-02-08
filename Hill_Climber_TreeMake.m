clc;
clear all;

format long
tic

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
%%%%%%%%%%%%%%%%%%CreatInitial%%%%%%%%%%%%%%%%%
matrix(1)=fuhao(1);

matrix(2)=fuhao(2); 

matrix(3)=fuhao(3); 

matrix(4)=fuhao(4); 
matrix(5)=fuhao(5); 
matrix(6)=fuhao(6); 
matrix(7)=99999;

for ii=[ceng4 ceng5]
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
xhowmany=randperm(20,1);
for jj=1:xhowmany
    xpick=randperm(rowmat-8,1)+8;
    matrix(xpick)=99999;
end

savefilename=['.\generate1st.mat'];
save(savefilename)
