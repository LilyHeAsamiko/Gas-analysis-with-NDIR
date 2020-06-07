% Example, how to solve a function numerically.

%Create vector x
x=1:0.01:10;

%Function y(x) that cannot be solved analytically
y=x.^2.*log((sin(x).^2+1).*x)+1./x+log(x)+1e-4.*x.^5;

%Plot how the function looks like
figure(1);
plot(x,y);

%Solve x at the point y=200

YY=200;
XX=interp1(y,x,YY);

