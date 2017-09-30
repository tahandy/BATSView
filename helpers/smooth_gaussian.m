function ynew = smooth_gaussian(x,yold,fwhm)
    
    if(numel(x)~=numel(yold))
        error('[smooth_gaussian] Inputs are of unequal size!');
    end
    
    if(sum(size(yold)>1)>1)
        error('[smooth_gaussian] smooth_gaussian only operates on 1D data');
    end
    
    if(fwhm<=0.e0)
        ynew = yold;
        return;
    end
    
    ynew  = zeros(size(yold));
    sigma = fwhm/2.355;
    for i=1:numel(yold)
        la = abs(x-x(i))>2*fwhm;
        filt = exp(-(x(la)-x(i)).^2/(2*sigma).^2);
        fsum = sum(filt);
        ynew(i) = sum(yold(la).*filt)/fsum;
    end
    
end
    