    x=-4:5;
    y=1+2*(sin(2*pi*0.1*x+2)+0.3*randn(size(x)));%Sine + noise
    tic;
    [SineP]=sineFit(x,y)
    toc
    figure;
    xx=x(1):(x(2)-x(1))/10:x(end);%better resolution
    plot(x,y,xx,SineP(1)+SineP(2)*sin(2*pi*SineP(3)*xx+SineP(4)));
%uncomment following line if you want to save y=f(x) and run it sineFitDemo
% save('xy.mat','x','y');