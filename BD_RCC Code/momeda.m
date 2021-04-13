function [MKurt f y] = momeda(x,filterSize,window,period,plotMode)
    % MULTIPOINT OPTIMAL MINUMUM ENTROPY DECONVOLUTION ADJUSTED
    %       code by Geoff McDonald (glmcdona@gmail.com), 2015
    %       Used in my research with reference to paper:
    %
    % momeda(x,filterSize,window,period,plotMode)
    %  Multipoint Optimal Minimum Entropy Deconvolution (MOMEDA) computation algorithm. This proposed
    %  method solves the optmial solution for deconvolving a periodic train of impulses from a signal.
    %  It is best-suited in application to rotating machine faults from vibration signals, to deconvolve
    %  the impulse-like vibration associated with many gear and bear faults.
    %
    %  This method is derived in the Algorithm Reference section.
    %
    %
    % Algorithm Reference:
    %    Preparing to publish.
    %    G.L. McDonald, <others>, Multipoint Optimal Minimum Entropy Deconvolution and Convolution
    %    Fix: Application to Vibration Fault Detection, unpublished
    %
    % Inputs:
    %    x: 
    %       Signal to generate apply MOMEDA on. Generally this should be around the range
    %       of 1000 to 10,000 samples covering at least 5 rotations of the elements in the machine.
    % 
    %    filterSize:
    %       This is the length of the finite inpulse filter filter to 
    %       design. Generally a number on the order of 500 or 1000 is good, but may
    %       depend on the dataset length.
    % 
    %    window:
    %       This is the window that be convolved with the impulse train target. Generally, a
    %       rectangular window works well, eg [1 1 1 1 1].
    % 
    %    period:
    %       This is the periods to test as the spectrum x-axis. It should be a decimal range, like:
    %           range = 5:0.1:300;
    % 
    %    plotMode:
    %       If this value is > 0, plots will be generated of the iterative
    %       performance and of the resulting signal.
    %
    % Outputs:
    %    MKurt:
    %       The Multipoint Kurtosis of the Deconvolution result.
    %
    %    f:
    %       Optimal FIR filter designed.
    %
    %    y:
    %       Filtered output signal.
    %
    % Example:
    %
    % % Simple vibration fault model
    % close all
    % n = 0:4999;
    % h = [-0.05 0.1 -0.4 -0.8 1 -0.8 -0.4 0.1 -0.05];
    % faultn = 0.05*(mod(n,50)==0);
    % fault = filter(h,1,faultn);
    % noise = wgn(1,length(n),-25);
    % x = sin(2*pi*n/30) + 0.2*sin(2*pi*n/60) + 0.1*sin(2*pi*n/15) + fault;
    % xn = x + noise;
    % 
    % % No window. A 5-sample recangular window would be ones(5,1).
    % window = ones(1,1);
    % 
    % % 1000-sample FIR filters will be designed
    % L = 1000;
    % 
    % % Recover the fault signal of period 50
    % [MKurt f y] = momeda(xn,L,window,50,1);
    %
    %
    
    % Assign default values for inputs
    if( isempty(filterSize) )
        filterSize = 300;
    end
    if( isempty(plotMode) )
        plotMode = 0;
    end
    if( isempty(window) )
        window = ones(1,1);
    end
    
    if( sum( size(x) > 1 ) > 1 )
        error('MOMEDA:InvalidInput', 'Input signal x must be 1d.')
    elseif(  sum(size(plotMode) > 1) ~= 0 )
        error('MOMEDA:InvalidInput', 'Input argument plotMode must be a scalar.')
    elseif( sum(size(filterSize) > 1) ~= 0 || filterSize <= 0 || mod(filterSize, 1) ~= 0 )
        error('MOMEDA:InvalidInput', 'Input argument filterSize must be a positive integer scalar.')
    elseif( sum(size(window) > 1) > 1 )
        error('MOMEDA:InvalidInput', 'Input argument window must be 1d.')
    elseif( period <= length(window) )
        error('MOMEDA:InvalidInput', 'Period should be larger than the length of the window.')
    elseif( filterSize >= length(x) )
        error('MOMEDA:InvalidInput', 'Input argument filterSize must be smaller than the length of input signal x.')
    end
    
    L = filterSize;
    x = x(:); % A column vector
    
    %%% Calculte X0 matrix
    N = length(x);
    X0 = zeros(L,N);
    
    for( l =1:L )
        if( l == 1 )
            X0(l,1:N) = x(1:N);
        else
            X0(l,2:end) = X0(l-1, 1:end-1);
        end
    end
    
                        % "valid" region only
    X0 = X0(:,L:N-1);   % y = f*x where only valid x is used
                        % y = Xm0'*x to get valid output signal
    
    autocorr = X0*X0';
    autocorr_inv = pinv(autocorr);
    
    % Built the impulse train vector separated the by periods
    t = zeros(N-L,1);
    points{1} = 1:period:(size(X0,2)-1);
    points{1} = round(points{1});
    t(points{1},1) = 1;
    
    % Apply the windowing function to the target vectors
    t = filter(window, 1, t);
    
    % Calculate the spectrum of optimal filters
    f = autocorr_inv * X0 * t;

    % Calculate the spectrum of outputs
    y = X0'*f;
    
    % Calculate the spectrum of PKurt values for each output
    MKurt = mkurt(y,t);
    
    % Plot the result
    if( plotMode > 0 )
        figure;
        subplot(3,1,1)
        plot(x)
        title('Input signal');
        xlabel('Sample number');
        
        subplot(3,1,2)
        plot(y)
        title('Output signal');
        xlabel('Sample number');
        
        subplot(3,1,3)
        stem(f)
        title('Designed optimal FIR filter');
        xlabel('Sample number');
    end
end

function [result] = mkurt(x,target)
    % This function simply calculates the summed kurtosis of the input
    % signal, x, according to the target vector positions.
    result = zeros(size(x,2),1);
    for i = 1:size(x,2)
        result(i) = ( (target(:,i).^4)'*(x(:,i).^4) )/(sum(x(:,i).^2)^2) * sum(abs(target(:,i)));
    end
end
