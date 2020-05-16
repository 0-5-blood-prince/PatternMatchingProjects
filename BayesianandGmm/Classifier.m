close all;
clear all;
%type trian.txt
D = readmatrix("NLStrian.txt");
Dd=readmatrix("NLSdev.txt");
X = D(:,1:1);
Y = D(:,2:2);
P = D(:,1:2);
C = D(:,3:3);
%plot(X(1:350,:),Y(1:350,:),".")hold on
%plot(X(351:700,:),Y(351:700,:),"o")hold on
%plot(X(701:1050,:),Y(701:1050,:),"*")hold on
sx1=0;sy1=0;n1=0;
sx2=0;sy2=0;n2=0;
sx3=0;sy3=0;n3=0;
sxt=0;syt=0;
for i = 1:1050
    sxt=sxt+X(i,1);
    syt=syt+Y(i,1);
    if(C(i,1)==1)
        sx1= sx1+X(i,1);
        sy1= sy1+Y(i,1);
        n1=n1+1;
    elseif(C(i,1)==2)
         sx2=sx2+X(i,1);
         sy2=sy2+Y(i,1);
         n2=n2+1;
    else
        sx3=sx3+X(i,1);
        sy3=sy3+Y(i,1);
        n3=n3+1;
    end  
end
Me1(1,1)=sx1/n1;Me2(1,1)=sx2/n2;Me3(1,1)=sx3/n3;
Me1(1,2)=sy1/n1;Me2(1,2)=sy2/n2;Me3(1,2)=sy3/n3;
Met(1,1)=sxt/1050;Met(1,2)=syt/1050;
Mett=(sxt+syt)/2010;
Vart=0;
Cov1=zeros(2,2);Cov2=zeros(2,2);Cov3=zeros(2,2);Covt=zeros(2,2);
st=0;
stt=0;
for i=1:2
    for j=1:2
        s=0;
        st=0;
        for k=1:1050
            if(i==j)
                stt=stt+(P(k,i)-Mett)*(P(k,i)-Mett);
            end
            st=st+(P(k,i)-Met(1,i))*(P(k,j)-Met(1,j));
            if (C(k,1)==1)
                s=s+(P(k,i)-Me1(1,i))*(P(k,j)-Me1(1,j));
            end
        end
        Cov1(i,j) = s/350;
        Covt(i,j) = st/1050;
    end
end
Vart=stt/2010;
for i=1:2
    for j=1:2
        s=0;
        for k=1:1050
            if (C(k,1)==2)
                s=s+(P(k,i)-Me2(1,i))*(P(k,j)-Me2(1,j));
            end
        end
        Cov2(i,j) = s/350;
    end
end
for i=1:2
    for j=1:2
        s=0;
        for k=1:1050
            if (C(k,1)==3)
                s=s+(P(k,i)-Me3(1,i))*(P(k,j)-Me3(1,j));
            end
        end
        Cov3(i,j) = s/350;
    end
end
%%%%
%Covariance Modification
%{
Cov1=Covt;
Cov2=Covt;
Cov3=Covt;
%}
%{
 Cov1=Vart*eye(2);
 Cov2=Vart*eye(2);
Cov3=Vart*eye(2);
%}
%{
 Cov1=Covt;
Cov2=Covt;
Cov3=Covt;
Cov1(1,2)=0;Cov1(2,1)=0;
Cov2(1,2)=0;Cov2(2,1)=0;
Cov3(1,2)=0;Cov3(2,1)=0;
%}
%{
Cov1(1,2)=0;Cov1(2,1)=0;
Cov2(1,2)=0;Cov2(2,1)=0;
Cov3(1,2)=0;Cov3(2,1)=0;
%}
%{
plotgaussian(Me1,Cov1,1);
hold on
plotgaussian(Me2,Cov2,2);
hold on
plotgaussian(Me3,Cov3,3);

hold off
%}

%case1

g1 = discriminant(Cov1,Me1');
%fimplicit(g1)
g2 = discriminant(Cov2,Me2');
g3 = discriminant(Cov3,Me3');

%{
dispdiscriminant(Cov1,Me1',Cov2,Me2');
hold on
dispdiscriminant(Cov2,Me2',Cov3,Me3');
hold on
dispdiscriminant(Cov1,Me1',Cov3,Me3');
hold off
%}
plotconfusionm(g1,g2,g3,Dd);

%[To,Fo]=plotroc(g1,g2,g3,Dd);
function []=plotconfusionm(k1,k2,k3,Ds)
    T1=0;F12=0;F13=0;
    T2=0;F21=0;F23=0;
    T3=0;F31=0;F32=0;
    Te=zeros(300,3);
    T=zeros(300,3);
   
    for i=1:300
        d=[Ds(i,1);Ds(i,2)];
        c1=k1(d);
        c2=k2(d);
        c3=k3(d);
       if(Ds(i,3)==1) 
           T(i,1)=1;
       end
       if(Ds(i,3)==2) 
           T(i,2)=1;
       end
       if(Ds(i,3)==3) 
           T(i,3)=1;
       end
        if((c1>c2) & (c1>c3))
            T1=T1+1;
            Te(i,1)=1;
        elseif((c2>c1) & (c2>c3))
            T2=T2+1;
            Te(i,2)=1;
        elseif((c3>c1) & (c3>c2))
            T3=T3+1;
            Te(i,3)=1;
        end
    end
    %{
targets = zeros(3,300);
outputs = zeros(3,300);
targetsIdx = sub2ind(size(targets), T, 1:300);
outputsIdx = sub2ind(size(outputs), Te, 1:300);
targets(targetsIdx) = 1;
outputs(outputsIdx) = 1;
    figure;
    plotconfusion(targets,outputs); 
  %}
    %disp(T-Te);
    figure;
   %plotconfusion(T,Te);
    cm=confusionchart(T,Te);
end
function [TPR1,FPR1]=plotroc(k1,k2,k3,Ds)
      TPR1=[];TPR2=[];
      FPR1=[];
      T=0;
      for t=0:1:100
          thres=(100-t)/100;
          TP1=0;TP2=0;TP3=0;
          
          FN1=0;FN2=0;FN3=0;
       n1=0;n2=0;n3=0;
      for i=1:300 
        d=[Ds(i,1);Ds(i,2)];
        c1=k1(d);
        c2=k2(d);
        c3=k3(d);T=Ds(i,3);
         l=0;k=0;m=0;
            if((c1/(c1+c2+c3)>thres)) 
               l=1;n1=n1+1;
            end
        
            if((c2/(c1+c2+c3)>thres))
               k=1;n2=n2+1;
            end
    
            if((c3/(c1+c2+c3)>thres))
               m=1;n3=n3+1;
            end
      
          if(l==1 & T==1) TP1=TP1+1;
          elseif(T==1)    FN1=FN1+1;
          end
      
       
          if(k==1 & T==2) TP2=TP2+1;
          elseif(T==2)     FN2=FN2+1;
          end
       
      
          if(m==1 & T==3) TP3=TP3+1;
          elseif(T==3)     FN3=FN3+1;
          end
      
      end
       FP1=n1-TP1;FP2=n2-TP2;FP3=n3-TP3;
       TN1=200-FP1;TN2=200-FP2;TN3=200-FP3;
       TPR1 =[TPR1 ,TP1/(TP1+FN1)];
       FPR1 = [FPR1,FP1/(FP1+TN1)];
      end
      figure;
    
      plot(TPR1,FPR1);
end
function [K] = gaussianf(I,J,mean,covariancem)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
d = (-1/2)*[I J]*inv(covariancem)*[I J]';
modulo = det(covariancem);
K = (1/(2*pi*sqrt(modulo)))*exp(d);
end
function [] = plotgaussian(mean,covariancem,c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Zp = zeros(100,100);
 x = linspace(-1000,100,100) ;
 y = linspace(-2500,100,100) ;
[Xp,Yp] = meshgrid(x,y);
 for i=1:100
     for j=1:100
         Zp(i,j)=gaussianf(Xp(i,j)+mean(1,1),Yp(i,j)+mean(1,2),mean,covariancem);
     end
 end
 
%surfc(Xp,Yp,Zp);
zlim([-0.00001,0.00001]);

 contour(Xp,Yp,Zp);
 hold on;
 [U,V]=eig(covariancem);
 syms x1 x2;
 p = [x1+mean(1,1);x2+mean(1,2)];
 e1=U(1,:)*p;
 e2=U(2,:)*p;
disp(vpa(e1,3));
disp(vpa(e2,3));
range=mean(1,1);
range2=mean(1,2);
ezplot(e1,[-1000 100 -2500 100]);
hold on
ezplot (e2,[-1000 100 -2500 100]);
hold on


end
function [g]=discriminant(cov,mean)
W = (-1/2)*inv(cov);
w1 = inv(cov)*mean;
wo=-(1/2)*(mean)'*inv(cov)*(mean)-(0.5)*log(det(cov));
g =@(x) x'*W*x+w1'*x+wo;
end
function [] =dispdiscriminant(cov1,mean1,cov2,mean2)
W1 = (-1/2)*inv(cov1);
w11 = inv(cov1)*mean1;
wo1=-(1/2)*(mean1)'*inv(cov1)*(mean1)-(0.5)*log(det(cov1));
syms x y;
p = [x;y];
axis([-30,40,-30,40]);
g1=p'*W1*p+w11'*p+wo1;
disp('g1')
disp(vpa(g1,3));
W2 = (-1/2)*inv(cov2);
w12 = inv(cov2)*mean2;
wo2=-(1/2)*(mean2)'*inv(cov2)*(mean2)-(0.5)*log(det(cov2));
g2=p'*W2*p+w12'*p+wo2;
disp('g2');
disp(vpa(g2,3));
disp('g1-g2');
disp(vpa(g1-g2,3))
ezplot(g1-g2,[-30 40 -30 40]);
end