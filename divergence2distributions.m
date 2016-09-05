function [kd, jd] = divergence2distributions(P, Q)
%P = rand(1,10000); Q = rand(1,10000);
% N = 100; p = 0.5; P = randsrc(N,N,[0 1; (1-p) p]); Q = randsrc(N,N,[0 1; (1-p) p]);
%convert P and Q binary matrices into one dimensional vector
lp = length(P); lq= length(Q);
P = reshape(P, 1, lp*lp);
Q = reshape(Q, 1, lq*lq);
[hp, ep] = histcounts(P,2); [hq,eq] = histcounts(Q,2);
kd = KLDiv(hp,hq); jd = JSDiv(hp,hq);
%[JSD,Smax]=DJS_and_significanse(hp,hq,2);
fprintf('KLDiv= %s, JSDiv=%s, JSD=%s Smax=%s\n', num2str(kd), num2str(jd));%, num2str(JSD), Smax);
end


% for i=1:lp
%     if P(i) > 0 P(i) = 2;
%     elseif P(i) > 0.5 P(i) = 1;
%     else P(i)=0;
%     end
%       if Q(i) > 0.6 Q(i) = 2;
%       elseif Q(i) > 0.1  Q(i) = 1;
%       else Q(i)=0;
%       end
% end
% hp = hist(P,3); hq = hist(Q,3);
% kd = KLDiv(hp,hq); jd = JSDiv(hp,hq);  