function nL = compute_patch_each_layer(slip_model)
% compute the patches of each layer of each segment
% each row corresponds to each fault segment

   fault_id = slip_model(:,1);
   indx_layer = slip_model(:,3);
   
   Nf = max(fault_id);                      % number of fault segments
   Nl = max(indx_layer);                    % number of layers  
   Wp = zeros(1,Nf);                        % number of layers of each segment
   nL = zeros(Nf,Nl);                       % number of patches in each segment/each layer
   
   for ii = 1:Nf
       indx_this_segment = fault_id == ii;
       layer_this_segment = indx_layer(indx_this_segment);
       Wp(ii) = max(layer_this_segment);
   end
   
   for ii = 1:Nf
       indx_this_segment = fault_id == ii;
       layer_this_segment = indx_layer(indx_this_segment);
       for jj = 1:Wp(ii)
           tmp = find(layer_this_segment == jj);
           nL(ii,jj) = length(tmp);
       end
   end   
end