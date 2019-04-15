function [aic]= AIC(signal)

    Nn=length(signal);

    aic =zeros(Nn,1);
    
    % could be improved and vectorized 
    
    for k=1:length(aic)
       
        aic(k) =k*log(var(signal(1:k)))+(Nn-k-1)*log(var(signal(k+1:end)));
    end
    
end
