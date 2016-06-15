clear all;
close all;
% I=imread('lena_noisy.png');
I=imread('lena_noisy.jpg');
I = im2double(I(:,:,1));
maxpic=I(1,1);
minpic=I(1,1);
%values of I between 0 and 1
for i=1:length(I)
    for j=1:length(I)
        if(I(i,j)>maxpic)
            maxpic=I(i,j);
        end
        if(I(i,j)<minpic)
            minpic=I(i,j);
        end
    end
end

I=I/(maxpic-minpic);
Iold=I;
figure
imshow(Iold);

% steps=2000;
% 
% 
% %%%% for C=0.5*x^2, central difference
% %we'll use dx=dy=d=1
% 
% d=1;
% %CFL : dt<0.5*d^2
% dt=0.01*d^2;
% newI=I;
% % for C=0.5*x^2
% E=zeros(1,steps+1);
%     for p=1:length(I)
%         for q=1:length(I)
%     E(1)=E(1)+I(p,q)^2;
%         end
%         E(1)=0.5*E(1);
%     end
% 
% 
% for m=1:steps
%     % inside the picture
% for i=2:length(I)-1
%     for j=2:length(I)-1
%     newI(i,j)=I(i+d,j)-2*I(i,j)+I(i-d,j);
%     newI(i,j)=I(i,j+d)-2*I(i,j)+I(i,j-d);
%     newI(i,j)=newI(i,j)*dt/d+I(i,j);
%     end
% end
%     %edges
%     %i=1
%     for j=2:length(I)-1
%     newI(1,j)=I(1+d,j)-I(1,j);
%     newI(1,j)=I(1,j+d)-2*I(1,j)+I(1,j-d);
%     newI(1,j)=newI(1,j)*dt/d+I(1,j);
%     end
%     %i=length(I)
%     for j=2:length(I)-1
%     newI(length(I),j)=I(length(I),j)+I(length(I)-d,j);
%     newI(length(I),j)=I(length(I),j+d)-2*I(length(I),j)+I(length(I),j-d);
%     newI(length(I),j)=newI(length(I),j)*dt/d+I(length(I),j);
%     end
%     %j=1
%     for i=2:length(I)-1
%     newI(i,1)=I(i+d,1)-2*I(i,1)+I(i-d,1);
%     newI(i,1)=I(i,1+d)-I(i,1);
%     newI(i,1)=newI(i,1)*dt/d+I(i,1);
%     end
%     %j=length(I)
%     for i=2:length(I)-1
%     newI(i,length(I))=I(i+d,length(I))-2*I(i,length(I))+I(i-d,length(I));
%     newI(i,length(I))=-I(i,length(I))+I(i,length(I)-d);
%     newI(i,length(I))=newI(i,length(I))*dt/d+I(i,length(I));
%     end
%     
%     %computation of energy
%     for p=1:length(I)
%         for q=1:length(I)
%     E(m+1)=E(m+1)+I(p,q)^2;
%         end
%         E(m+1)=0.5*E(m+1);
%     end
%     
%     if m==1
%         I1=newI(100,:);
%     end
%     if m==200
%         I2=newI(100,:);
%     end
%     if m==300
%         I3=newI(100,:);
%     end
%      if m==2000
%         I4=newI(100,:);
%     end
%     
% I=newI;
% end
% 
% figure
% imshow(newI)
% figure
% plot([1:steps+1],E);
% xlabel(['Number of steps 0<x<' num2str(steps)]);
% ylabel('Energy');
% title('Perona-Malik diffusion  with c(x)=0.5*x^2')
% 
% 
% figure
% plot([1:215],I2,'b');
% hold on
% plot([1:215],I3,'g');
% hold on
% plot([1:215],I4,'r');
% xlabel(['y']);
% ylabel('I');
% title('I along the straight line x=100');

%for c(x)=x with hack, central difference
eps=0.00001;
%CFL condition
D=1;
Dt=sqrt(eps)*D^2*0.1;
In=I;
Ix=0; % ajouter le pas en divisions diff de 1
Iy=0;
Ixy=0;
stepz=400;
Eh=zeros(1,stepz);  
li=length(I);
lj=li;
lap=0;
denom=0;
for m=1:stepz
    Eh(1,m)=0;
for i=2:length(I)-1
    for j=2:length(I)-1
        Ix=0.5*(I(i+1,j)-I(i-1,j)); %central difference
        Iy=0.5*(I(i,j+1)-I(i,j-1));
        Ixy=0.25*(I(i+1,j+1)-I(i+1,j-1)-I(i-1,j+1)+I(i-1,j-1));
        In(i,j)=Ix^2*(I(i,j+1)-2*I(i,j)+I(i,j-1));
        In(i,j)=In(i,j)-2*Ix*Iy*Ixy;
        In(i,j)=In(i,j)+Iy^2*Ixy;
        denom=(Ix^2+Iy^2+eps)*sqrt((Ix^2+Iy^2+eps));
        lap=(I(i+1,j)-4*I(i,j)+I(i-1,j)+I(i,j+1)+I(i,j-1));
        In(i,j)=In(i,j)+eps*lap;
        In(i,j)=Dt*In(i,j)/denom + I(i,j);
        Eh(1,m)=Eh(1,m)+sqrt(Ix^2+Iy^2+eps);
    end
end
%     %i=1
%     for j=2:length(I)-1
%         Ix=0.5*(I(2,j)-I(1,j)); %central difference
%         Iy=0.5*(I(1,j+1)-I(1,j-1));
%         Ixy=0.25*(I(2,j+1)-I(2,j-1)-I(1,j+1)+I(1,j-1));
%         In(1,j)=Ix^2*(I(1,j+1)-2*I(1,j)+I(1,j-1));
%         In(1,j)=In(1,j)-2*Ix*Iy*Ixy;
%         In(1,j)=In(1,j)+Iy^2*Ixy;
%         In(1,j)=In(1,j)+eps*(I(2,j)-4*I(1,j)+I(1,j)+I(1,j+1)+I(1,j-1));
%         In(1,j)=Dt*In(1,j)/( (Ix^2+Iy^2+eps)*sqrt((Ix^2+Iy^2+eps)) ) + I(1,j);
%         Eh(1,m)=Eh(1,m)+sqrt(Ix^2+Iy^2+eps);
%     end
%     %i=length(I)
%     for j=2:length(I)-1
%         Ix=0.5*(I(li,j)-I(li-1,j)); %central difference
%         Iy=0.5*(I(li,j+1)-I(li,j-1));
%         Ixy=0.25*(I(li,j+1)-I(li,j-1)-I(li-1,j+1)+I(li-1,j-1));
%         In(li,j)=Ix^2*(I(li,j+1)-2*I(li,j)+I(li,j-1));
%         In(li,j)=In(li,j)-2*Ix*Iy*Ixy;
%         In(li,j)=In(li,j)+Iy^2*Ixy;
%         In(li,j)=In(li,j)+eps*(I(li,j)-4*I(li,j)+I(li-1,j)+I(li,j+1)+I(li,j-1));
%         In(li,j)=Dt*In(li,j)/( (Ix^2+Iy^2+eps)*sqrt((Ix^2+Iy^2+eps)) ) + I(li,j);
%         Eh(1,m)=Eh(1,m)+sqrt(Ix^2+Iy^2+eps);
%     end
%     %j=1
%     for i=2:length(I)-1
%         Ix=0.5*(I(i+1,1)-I(i-1,1)); %central difference
%         Iy=0.5*(I(i,2)-I(i,1));
%         Ixy=0.25*(I(i+1,2)-I(i+1,1)-I(i-1,2)+I(i-1,1));
%         In(i,1)=Ix^2*(I(i,2)-2*I(i,1)+I(i,1));
%         In(i,1)=In(i,1)-2*Ix*Iy*Ixy;
%         In(i,1)=In(i,1)+Iy^2*Ixy;
%         In(i,1)=In(i,1)+eps*(I(i+1,1)-4*I(i,1)+I(i-1,1)+I(i,2)+I(i,1));
%         In(i,1)=Dt*In(i,1)/( (Ix^2+Iy^2+eps)*sqrt((Ix^2+Iy^2+eps)) ) + I(i,1);
%         Eh(1,m)=Eh(1,m)+sqrt(Ix^2+Iy^2+eps);
%     end
%     %j=length(I)
%     for i=2:length(I)-1
%         Ix=0.5*(I(i+1,lj)-I(i-1,lj)); %central difference
%         Iy=0.5*(I(i,lj)-I(i,lj-1));
%         Ixy=0.25*(I(i+1,lj)-I(i+1,lj-1)-I(i-1,lj)+I(i-1,lj-1));
%         In(i,lj)=Ix^2*(I(i,lj)-2*I(i,lj)+I(i,lj-1));
%         In(i,lj)=In(i,lj)-2*Ix*Iy*Ixy;
%         In(i,lj)=In(i,lj)+Iy^2*Ixy;
%         In(i,lj)=In(i,lj)+eps*(I(i+1,lj)-4*I(i,lj)+I(i-1,lj)+I(i,lj)+I(i,lj-1));
%         In(i,lj)=Dt*In(i,j)/( (Ix^2+Iy^2+eps)*sqrt((Ix^2+Iy^2+eps)) ) + I(i,lj);
%         Eh(1,m)=Eh(1,m)+sqrt(Ix^2+Iy^2+eps);
% end
        if m==1
        I1=In(100,:);
        end
    I=In;
    
        if m==200
            I200=In(100,:);
        end
        
        if m==100
            I100=In(100,:)
        end
        
        
%     %normalization
%     maxpic=I(1,1);
%     minpic=I(1,1);
%     for i=1:length(I)
%     for j=1:length(I)
%         if(I(i,j)>maxpic)
%             maxpic=I(i,j);
%         end
%         if(I(i,j)<minpic)
%             minpic=I(i,j);
%         end
%     end
% end
% 
% I=I/(maxpic-minpic);
end
figure
imshow(In)
figure
plot([1:stepz],Eh);
figure

Ifin=In(100,:)
Iinit=Iold(100,:);
plot([1:215],Iinit,'r');
hold on
plot([1:215],Ifin,'b');
hold on
plot([1:215],I100,'g');
hold on
plot([1:215],I200,'k');
