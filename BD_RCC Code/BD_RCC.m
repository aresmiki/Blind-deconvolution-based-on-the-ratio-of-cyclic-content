function [optW,rec,funv]=BD_RCC(x,np,fault,Fs)
% Blind deconvolution based on RCC
    %  code by Liu He (aremiki@163.com), July 2020
    %  Used in my PhD research at the University of SouthWest Jiaotong University.
    %
    %    Algorithm Reference:
    %    [1] XXXXXXXXXXX.
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
    %       The fault period 
    %
    %    Fs: 
    %       The sampling frequency.
    %
    % Outputs:
    %    optW:
    %         The final 1d filter in finite impulse response format.
    %
    %    rec:
    %       The input signal(s) x, filtered by the resulting our method's filter.
    %
    %    funv:
    %       the final RCC vaule 
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

% x-信号
%% 构造Hankel 矩阵
x=x-mean(x);
N=length(x);
NP=floor(N-np+1);
data=zeros(NP,np);
for i=1:NP
    data(i,:)=x(i:(i+np-1));
end

%% 根据Fault计算特征谐波对应的向量矩阵 注意该Fault需要给的是包络谱中第一个特征频率

% tzf=floor(fault./(Fs./NP));
def=(Fs/NP);
faultX=[];
for i=1:3 % 前3阶谐波足够
    tzf=round(i*fault./(Fs./NP))+1;
    faultX=[faultX,[tzf-ceil(3/def):tzf+ceil(3/def)]]; %  循环平稳成分 窄带
end
%% 初始化权重
optW = lpc(x,np-1);
optW = optW(:);
optW=optW./norm(optW);

%%
Psi=ones(length(faultX),NP);
for i=1:length(faultX)
    Psi(i,:)=exp(-j*2*pi*(0:1:NP-1)*(faultX(i)-1)/NP);
end
[HM]=get_hilbfir_M(NP);  % Hilbert变换矩阵 Hilbert Matrix

% 优化
% 注意这里最小化的是-RCC，等价最大化RCC
% Note that the minimized -RCC is equivalent to the maximized RCC
[optW,funv,~,Info]  = minFunc(@RCC,optW(:), ... 
                   struct('MaxIter', 100,'TolX',1e-4,'Display',0),data,Psi,HM);

optW=optW./norm(optW);        
rec= data*optW;
end
