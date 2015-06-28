function csd = compute_csd_de(degvec)

% 0. scael to [-pi, pi]
radvec = degtorad(degvec);
% 1. deg to rad

% 2. compute the length
x = mean(sin(radvec));
y = mean(cos(radvec));
[~,R] = cart2pol(x,y);
    
% 3. csd
csd = radtodeg(sqrt(-2*log(R)));

function rad = degtorad(deg)
rad = deg*pi/180;

function deg = radtodeg(rad)
deg = rad*180/pi;

