%%
df = 5;
x = trnd(df, 400, 1);
x = x/std(x);
fit = fitdist(x, 'tlocationscale')
var(x)
%%
lambda = sampling_mixture_scale_parameter(x,1,fit.nu);
%%
dfs = ones(1000,1)*10;
for i=2:1000
    if mod(i,100) == 0
        disp(i)
    end
    lambda = sampling_mixture_scale_parameter(x,1,dfs(i-1));
    dfs(i) = sampling_degree_of_freedom(lambda,dfs(i-1),5,1);
end
plot(dfs)