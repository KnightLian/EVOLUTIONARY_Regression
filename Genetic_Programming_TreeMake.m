clc;
clear all;

format long
tic

evasize=30;%Must Even
rows=63;
ii=1;
evamatrix=zeros(evasize,rows+1);%+1 is distance

while ii<=evasize
    random63=new63mat();
    randtrans=random63.';
    evamatrix(ii,1:rows)=evamatrix(ii,1:rows) + randtrans;  
    ii=ii+1;   
end

savefilename=['.\generate30.mat'];
save(savefilename)

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
    randmat(1)=fuhao(randperm(6,1));
    routput63=randmat;    
end
