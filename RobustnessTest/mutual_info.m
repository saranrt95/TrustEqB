function MI=mutual_info(x,y)
    
    % joint probability
    pxy_obj=histogram2(x,y,'Normalization','probability');
    
    % per non plottare istogrammi
    % Obtain the handle of the figure that contains the histogram
    handleOfHistogramFigure = ancestor(pxy_obj, 'figure');
    % Make the figure window visible in case it was invisible before
    handleOfHistogramFigure.Visible  = 'off';
    
    pxy=nonzeros(pxy_obj.Values);
    %joint entropy (usando log2 l'entropia è in unità bits)
    Hxy=-dot(pxy,log2(pxy));
    
    % marginal probabilities
    Px_obj=histogram(x,'Normalization','probability');
    
    % per non plottare istogrammi
    % Obtain the handle of the figure that contains the histogram
    handleOfHistogramFigure = ancestor(Px_obj, 'figure');
    % Make the figure window visible in case it was invisible before
    handleOfHistogramFigure.Visible  = 'off';
    
    Px=nonzeros(Px_obj.Values);
    
    Py_obj=histogram(y,'Normalization','probability');
    
    % per non plottare istogrammi
    % Obtain the handle of the figure that contains the histogram
    handleOfHistogramFigure = ancestor(Py_obj, 'figure');
    % Make the figure window visible in case it was invisible before
    handleOfHistogramFigure.Visible  = 'off';
    
    Py=nonzeros(Py_obj.Values);
    % marginal entropies
    Hx=-dot(Px,log2(Px));
    Hy=-dot(Py,log2(Py));
    MI=Hx+Hy-Hxy;
    
end

% provare con muInf += pab2d[j][i] * log(pab2d[j][i]/p1[i]/p2[j]);