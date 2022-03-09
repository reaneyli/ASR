function [xx] = update_marker(x_sub,x_or,w,d)
  [n,dim] = size(x_or); 
  [n_sub,dim] = size(x_sub); 
  xx=x_sub;
  for i=1:n_sub
      x = [ones(1,1), x_sub(i,:)];
      K= zeros (1,n);  
      for j=1:n
          K(1,j)=sum((x_sub(i,:)-x_or(j,:)).^2,2);       
      end    
      K = - sqrt(K);
%        mask = K < 1e-10; % to avoid singularity.    
%       %K = - sqrt(K).* (K>1e-10);                % For Face3D
%       K = 0.5 * K .* log(K + mask) .* (K>1e-10); % For 2D Demo cases
      xx1 = x*d + K*w;
      xx(i,:) = xx1(:,2:dim+1);
  end
end