function [c,D,count] = kmc(x,cn,plot_ch)
% KMC    [c,D,count]=kmc(x,cn)  K-means clustering
%    x:        N-by-L, contains L number of L-length data vectors to be clustered 
%    cn:       1-by-1, preselected number of clusters
%    c:        N-by-cn, cluster centers
%    D:        N-by-cn, covariance matrix of each cluster
%    count:    1-by-cn, count(k) is the number of data in the k-th cluster
%    code:     1-by-L, optional, the index of class to which each data vector belongs, 
%              all codes constitue a codebook in HMM
%
% Copyright (C) 2000-2012 by Xuejun Liao @ Duke University
% Send comments to: xuejun.liao@gmail.com, xjliao@duke.edu, (919)-660-5548
% Written by Xuejun Liao,  June 15, 2000


[N,L]=size(x);

% Select the initial centers by maximizing the between-center distances 
index = 1; %ceil(cn*rand);  
c         = x(:,index); 
for i = 2:cn
    dist = zeros(1,L);
   for j=1:i-1
      temp_x = x; 
      temp_x(:,index) = c(:,j)*ones(1,i-1);
      dist = dist + sum((temp_x-c(:,j)*ones(1,L)).^2,1);
   end
   [aa index(i)] = max(dist);
   c(:,i) = x(:,index(i));
end

% K-means clustering 
t=1;  err=1;   last_class_index=zeros(1,L);
while err %>10^(-9)  
   % assign each input data vector to a cluster by minimum-distance criterion
   class_index=0;  d_tol=[];
   for ci=1:cn
      xdif = x-repmat(c(:,ci),[1,L]);
      dist(ci,:) = sum((x-repmat(c(:,ci),[1,L]) ).^2,1);      
   end
   [aa class_index]=min(dist,[],1);       
   last_c = c;
   idx = unique(class_index);
   if length(idx)<size(c,2)
       c = c(:,unique(class_index));
       cn = length(idx);
       class_index=0;  d_tol=[];
       for ci=1:cn
           xdif = x-repmat(c(:,ci),[1,L]);
           dist(ci,:) = sum((x-repmat(c(:,ci),[1,L]) ).^2,1);  
       end
       [aa class_index]=min(dist,[],1);
       last_c = c;   
   end
   
   % adjust the centers 
   for i=1:cn
      c(:,i)=mean(x(:,find(class_index==i)),2);
   end
   
   % examine the termination condition   
   %err=sum(sum(((c-last_c).*(c-last_c)).^2))/sum(sum(((c+last_c).*(c+last_c)).^2));
   err=length(find(last_class_index~=class_index));   
   if strcmp(plot_ch,'on')
       fprintf(1,'K-means clustering --- %d-th iteration, error=%f\n',t,err);
   end
   t = t+1; 
   last_class_index=class_index;  
end

code   = class_index;
D      = cell(1,cn); 
count = [];  
for i = 1:cn
   index = find(class_index==i);
   len = length(index);
   if len>0
       xdif = (x(:,index) - repmat(c(:,i),[1,len]));
       D{i} = xdif*xdif'/length(index);
       D{i} = D{i} + eye(N)*1e-5*max(eig(D{i}));
   end         
   count(i) = length(index);
end
idx = find(count>=2); %floor(L/cn/2));
c = c(:,idx);
D = D(idx);
count = count(idx);
count = count/sum(count);