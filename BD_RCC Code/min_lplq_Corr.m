%% 最小化L2/L1
function [optW,rec,funv]=min_lplq_Corr(x,np,fault,p,q)
%%
    % Blind deconvolution based on CG-Lp/Lq
    %  code by Liu He (aremiki@163.com), July 2020
    %  Used in my PhD research at the University of SouthWest Jiaotong University.
    %
    %    Algorithm Reference:
    %    [1] L. He, Y. Li, Y. Liu and J. Lin, "Minimum Correlated Generalized Lp/Lq Deconvolution 
    %         for Recovering Repetitive Impacts From a Vibration Mixture," in IEEE Sensors Journal, 
    %         vol. 21, no. 2, pp. 2043-2054, 15 Jan.15, 2021, doi: 10.1109/JSEN.2020.3021213.
    %
    % Inputs:
    %    x: 
    %       Signal to perform blind deconvolution on. 
    % 
    %    np:
    %       This is the length of the finite inpulse filter filter to 
    %       design. Investigate the performance difference using
    %       different values.
    %    fault:
    %       The fault period (in sample)
    %
    %    p and q: (OPTIONAL)
    %       The Sparsity criterion CG-Lp/Lq, CG-L1/L2 is the default value.
    %
    % Outputs:
    %    optW:
    %         The final 1d filter in finite impulse response format.
    %
    %    rec:
    %       The input signal(s) x, filtered by the resulting our method's filter.
    %
    %    funv:
    %       the final CG-Lp/Lq 
    % 
    % Note 1:
    %    When using this code, add the minFunc file to the MATLAB environment
    %    Execute the code "addpath('./minFunc');"in this path.
    %
    % Note 2:
    %   The solution is not guaranteed to be the optimal solution 
    %   to our minimization problem, the solution is just a local 
    %   minimum and therefore a good pick. 
    %
    % Note 3:
    %     If you want to get the global optimal solution, please refer to 
    %     the optimization techniques in the following literature
    %    [1] Liu He, et.al, Optimized minimum generalized Lp/Lq deconvolution for recovering repetitive 
    %    impacts from a vibration mixture, Measurement, vol. 168，
    %     https://doi.org/10.1016/j.measurement.2020.108329 
    
if nargin <= 3
    p=1;q=2;
end
% x-信号
% np-神经元个数
%% 构造Hankel 矩阵
N=length(x);
data=[];
for i=1:floor(N-np+1)
    data(i,:)=x(i:(i+np-1));
end
%% 初始化权重
optW = zeros(size(data, 2),1); %使用高峰度的单点函数
optW(2)=1;
optW=optW./norm(optW);
%%
[optW,~,~,funv] = minFunc(@Corr_obj_grad, optW(:), ...
                   struct('MaxIter', 400,'Display','final'), data,p,q,fault);
optW=optW./norm(optW);        
rec= data*optW;
end