function [grid,data,fit,vals] = pdffit_modcau(a,guess_a,guess_b,guess_C)
% Fits a modified Cauchy distribution to the data a

[data,edges] = histcounts(a,1e2,'Normalization','pdf');
grid = 0.5*(edges(1:end-1)+edges(2:end));

p = @(a,x) a(3)*(a(1)^2+x.^2).^(-a(2)/2);

options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',2000);
vals = lsqcurvefit(p,[guess_a,guess_b,guess_C],grid,data,[0,1,0],[],options);

fit = p(vals,grid);

end