function [x,w] = Gaus_int(n)

%%Works till n == 5

if n == 2
    x = [-0.5773502691896257 0.5773502691896257]';
    w = [1 1]';
elseif n==3
    x = [0 -0.7745966692414834 0.7745966692414834]';
    w = [0.8888888888888888 0.5555555555555556 0.5555555555555556]';
elseif n==4
    x = [-0.3399810435848563 0.3399810435848563 -0.8611363115940526 0.8611363115940526]';
    w = [0.6521451548625461 0.6521451548625461 0.3478548451374538 0.3478548451374538]';
elseif n==5
    x = [0 -0.5384693101056831 0.5384693101056831 -0.9061798459386640 0.9061798459386640]';
    w = [0.5688888888888889 0.4786286704993665 0.4786286704993665 0.2369268850561891 0.2369268850561891]';
end
end
