
function [tempo,mrp,mrp_dot,mrp_2dot] = bspline_interpolation_simple(u,tf,mrp1,mrp2,mrp3,N0,N1,N2,t,k,m)


nn = m-1;
mrp = zeros(3,size(N0,1));
mrp_dot = zeros(3,size(N0,1));
mrp_2dot = zeros(3,size(N0,1));
tempo = u*tf;

%%% COSTRUZIONE DELLE CURVE CON LE B-SPLINES
mrp(1,:) = N0*mrp1';
mrp(2,:) = N0*mrp2';
mrp(3,:) = N0*mrp3';

%%% Derivata prima
qy1 = zeros(nn,1);
qy2 = zeros(nn,1);
qy3 = zeros(nn,1);
tt = t(2:end-1);

for ii = 1:nn
    qy1(ii) = (k-1)*(mrp1(ii+1) - mrp1(ii))/(tt(k+ii-1)-tt(ii));
    qy2(ii) = (k-1)*(mrp2(ii+1) - mrp2(ii))/(tt(k+ii-1)-tt(ii));
    qy3(ii) = (k-1)*(mrp3(ii+1) - mrp3(ii))/(tt(k+ii-1)-tt(ii));    
end
f1_dot = N1*qy1;
f2_dot = N1*qy2;
f3_dot = N1*qy3;

mrp_dot(1,:) = f1_dot/tf;
mrp_dot(2,:) = f2_dot/tf;
mrp_dot(3,:) = f3_dot/tf;

%%% Derivata seconda

qyy1 = zeros(nn-1,1);
qyy2 = zeros(nn-1,1);
qyy3 = zeros(nn-1,1);
ttt = t(3:end-2);

for ii = 1:nn-1
    qyy1(ii) = (k-2)*(qy1(ii+1) - qy1(ii))/(ttt(k+ii-2)-ttt(ii));
    qyy2(ii) = (k-2)*(qy2(ii+1) - qy2(ii))/(ttt(k+ii-2)-ttt(ii));
    qyy3(ii) = (k-2)*(qy3(ii+1) - qy3(ii))/(ttt(k+ii-2)-ttt(ii));
end

f1_2dot = N2*qyy1/tf^2;
f2_2dot = N2*qyy2/tf^2;
f3_2dot = N2*qyy3/tf^2;


mrp_2dot(1,:) = f1_2dot;
mrp_2dot(2,:) = f2_2dot;
mrp_2dot(3,:) = f3_2dot;

end