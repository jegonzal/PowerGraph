function plot_powerlaw(degree)

counts = histc(degree, 1:max(degree));
nzind = counts > 0;
x = (1:max(degree))';
y = counts;
nzx = x(nzind);
nzy = y(nzind);
set(gcf, 'Position',[0,0,800,600]);
set(gca,'box','on');    
set(gca, 'FontSize', 18);
loglog(nzx, nzy, 'kx', 'LineWidth', 2);
xlabel('Degree');
ylabel('Count');



nzdegree = degree(degree > 0);
alpha = 1 + length(nzdegree) / sum(log(nzdegree / min(nzdegree)));
disp(['Alpha estimate: ', num2str(alpha)]);

hold on;
est = nzx.^(-alpha);
loglog(nzx, sum(nzy) / sum(est) * est, 'r','LineWidth', 2);
hold off;


% logx = [log(nzx), ones(size(nzx))];
% logy = log(nzy);
% w = regress(logy, logx);
% alpha = -w(1);
% disp(['Alpha regression estimate: ', num2str(alpha)]);
% 
% hold on;
% est = nzx.^(-alpha);
% loglog(nzx, sum(nzy) / sum(est) * est, 'b','LineWidth', 2);
% hold off;

end