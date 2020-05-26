function [grid,data,fit,vals] = pdffit_exp(a,guess_a,guess_b)
% Fits a distribution of the form p to the data a

[data,edges] = histcounts(a(:),5e2,'Normalization','pdf');
grid = 0.5*(edges(1:end-1)+edges(2:end));

p = @(a,x) a(2)*exp(a(2)*x-a(1)-exp(a(2)*x-a(1)));

options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',2000);
vals = lsqcurvefit(p,[guess_a,guess_b],grid,data,[],[],options);

fit = p(vals,grid);

end