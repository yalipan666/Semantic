function prob=sub_pca_dim(Data,Mask,N2)
% PCA_DIM   Automatic choice of dimensionality for PCA
%
% Laplace approximation and subsequernt simplifications to the
% posterior probability of the Data being generated from hidden var's
% of dimensionality k

% Methods described in
%   Thomas Minka:
%                  Automatic choice of dimensionality for PCA
%
%% Dimension of Data
% The Input Data MUST BE Channel x Freq x Time x Trial  
%                     or Channel x Time x Trial  
if ndims(Data)>2
    % If the dimension of Data is 3 or 4
    Data = Data(:,:);
end
% if (length(size(Data))==4)
%     Data=reshape(Data,size(Data,1)*size(Data,2)*size(Data,3), ...
%         size(Data,4))';
% end;
%%
[d,N]=size(Data);
E=0;
if d>1
    % if d>N
    %	 Data=Data';
    %	 [d,N]=size(Data)
    %      end;

    if (nargin==1)
        Mask=boolean(ones(1,N));
    end;

    if (length(size(Mask))==3)
        Mask=boolean(reshape(Mask,size(Mask,1)*size(Mask,2)*size(Mask,3),1)');
    end;


    Data=Data(:,Mask);
    [d,N]=size(Data);

    % Data=remmean(Data);
    % Data=remmean(Data')';
    % Data=Data./sqrt(var(Data')'*ones(1,N));
%     [E,D]=hyv_ica(Data,'only','pca','verbose','off');
%     lambda=diag(D)';
    [E,tmp1,lambda] = princomp(Data.','econ');
    lambda = lambda.';

    % [E,D]=eig(Data*Data'/size(Data,2));
    %lambda=fliplr(diag(D)');

    size(lambda);
    size((iFeta([0.05:0.01:4],size(lambda,2),N/2)));

    lambda = lambda./(iFeta([0.05:0.01:4],size(lambda,2),N/2));

else
    lambda=Data;
    N=N2;
    d=Mask;
end;



if prod(size(Mask))==1
    lambda=lambda(1:d);
end;
logLambda=log(lambda);
d=length(lambda);
k=[1:d];

m=d*k-0.5*k.*(k+1);

lgam=cumsum(gammaln((d-k+1)/2));

l_prob_U=-log(2).*k+lgam-cumsum(((d-k+1)/2)*log(pi));

tmp=fliplr(sum(lambda)-cumsum(lambda));
tmp(1)=min(abs(2*tmp(2)-tmp(3)),0.9*tmp(2));
tmp=fliplr(tmp);
tmp2=fliplr(sum(logLambda)-cumsum(logLambda));
tmp2(1)=tmp2(2);
tmp2=fliplr(tmp2);

tmp3=([tmp(1:d-1) tmp(d-1)])./([d-k(1:d-1) 1]);

l_nu=-N/2*(d-k).*log(tmp3);
l_nu(d)=0;

a=lambda;
b=tmp3;
% t1=cumsum(sum(log(triu((ones(d,1)*a)'-(ones(d,1)*a),1)),2));
% t2=cumsum(sum(log(triu(((b.^(-1))'*ones(1,d))'-(ones(d,1)*(a.^(-1)))',1)),2));
t1=triu((ones(d,1)*a)'-(ones(d,1)*a),1);
t2=triu(((b.^(-1))'*ones(1,d))'-(ones(d,1)*(a.^(-1)))',1);
t1(t1<=0)=1;
t2(t2<=0)=1;
t1=cumsum(sum(log(t1),2));
t2=cumsum(sum(log(t2),2));
l_Az=0.5*(t1+t1)';


l_lam=-(N/2)*cumsum(log(lambda));


l_prob_DK = (l_prob_U + l_nu + l_Az + l_lam + (m+k)/2*log(2*pi) - k/2*log(N));
l_BIC = (l_lam+l_nu-(m+k)/2*log(N));
l_RRN = (-N*k/2.*log(cumsum(lambda)./k)+l_nu);


% AIC and MDL

llhood = [1./(d-k(1:d-1))].*tmp2(1:d-1)-log([1./(d-k(1:d-1))].* tmp(1:d-1));

AIC = -2*N*(d-k(1:d-1)).*llhood + 2*(1+d.*k(1:d-1)+0.5*(k(1:d-1)-1));

MDL = -N*(d-k(1:d-1)).*llhood + 0.5*(1+d.*k(1:d-1)+0.5*(k(1:d-1)-1))*log(N);


% l_prob_DK=(l_prob_DK-min(l_prob_DK))./(max(l_prob_DK)-min(l_prob_DK));
% l_BIC=(l_BIC-min(l_BIC))./(max(l_BIC)-min(l_BIC));
% l_RRN=(l_RRN-min(l_RRN))./(max(l_RRN)-min(l_RRN));
% AIC=1-((AIC-min(AIC))./(max(AIC)-min(AIC)));
% MDL=1-((MDL-min(MDL))./(max(MDL)-min(MDL)));


prob.lap=0.95*l_prob_DK;
prob.bic=0.95*l_BIC;
prob.rrn=0.95*l_RRN;
prob.AIC=-0.95*AIC;
prob.MDL=-0.95*MDL;

prob.eig=lambda;
prob.leig=logLambda;

prob.lML=llhood;

prob.lAz=l_Az;
prob.llam=l_lam;
prob.lnu=l_nu;
prob.lpU=l_prob_U;

if size(E,1)>0
    prob.E=E;
end;

%% 
function [Result] = iFeta(eta,d1,d2);
   
  res = d1*(1-Feta(eta,d1/d2));
 
  Result = zeros(1,d1);
  
  for k= 1 : d1;
     Result(k) = eta(max(find(res>=k)));
  end;
  
%%   
  function [Result] = Feta(eta,y,sig2);
   
   if nargin<3
      sig2=1;
   end;
   
   bm=sig2*(1-sqrt(y))^2;
   bp=sig2*(1+sqrt(y))^2;
   
   Result=zeros(size(eta));

   teta=bm:1/1000:bp;
   feta=zeros(size(teta));
   feta=((2*pi*y*teta).^(-1)).*sqrt((teta-bm).*(bp-teta));
  
   tmp=(teta'*ones(1,size(eta,2)))./(ones(size(teta,2),1)*eta)<1;
   Result=sum((0.001*feta'*ones(1,size(eta,2))).*tmp);
  