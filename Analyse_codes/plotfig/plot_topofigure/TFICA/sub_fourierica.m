function [S_orig,A_orig,W_orig] = sub_fourierica(X, components, pcadim, complexmixing, removeoutliers, maxiter)
%SIMPLE CODE FOR DOING FOURIER-ICA
%as proposed in Hyvarinen, Ramkumar, Hari, Neuroimage 2010.
%Aapo Hyvarinen, Dec 2010

%Input data:
% X:                Input EEG/MEG data with dimension (channel x frequency x time x trial)
% components:       Number of independent components in decomposition
% pcadim:           Dimension after PCA on channels
% complexmixing:    Do we use complex-valued mixing (1) instead of real-valued (default 0)?
% removeoutliers:   Another option: Remove outliers? (default 1)
% maxiter:          Max number of iterations in estimation loop (default 1000)

%Output data:
% A_orig: spatial patterns
% W_orig: spatial filters
% S_orig: STFT's of the source signals

%% Extract parameters and re-organize data
[N_Channels,N_Freq,N_Time,N_Trials] = size(X);
X0 = single(X); clear X;
for n_freq=1:N_Freq; for n_chan=1:N_Channels;
	X(n_freq,:,n_chan) = reshape(X0(n_chan,n_freq,:,:),N_Time*N_Trials,1);
end; end;
clear X0;
N_MultipleTrialTime = N_Time*N_Trials;

%% Remove outliers, i.e. windows with really large norms
if removeoutliers
    fprintf('Outlier removal:');
    Xmat=reshape(X,[N_Freq,N_MultipleTrialTime*N_Channels]);
    lognorms=log(sum(abs(Xmat).^2));
    outlierthreshold=mean(lognorms)+3*std(lognorms); %what is good threshold?
    outlierindices=find(lognorms>outlierthreshold);
    fprintf(' removed %u windows\n',size(outlierindices,2));
    Xmat(:,outlierindices)=zeros(N_Freq,size(outlierindices,2));
    X=reshape(Xmat,[N_Freq,N_MultipleTrialTime,N_Channels]);
    clear Xmat
end

%% Do PCA and whitening wrt to channels (this is compulsory)
%  as typical in BSS algorithms

disp('Spatially whitening data with PCA dimension reduction');
%substract means for each channel
for k=1:N_Channels
    X(:,:,k)=X(:,:,k)-mean(mean(X(:,:,k)));
end
%"vectorize" data w.r.t. to last dimension (channels)
Xmat_c=(reshape(X,[N_Freq*N_MultipleTrialTime,N_Channels])).';
clear X
%store sample size after this transformation
N=size(Xmat_c,2);

%DO PCA AND WHITEN DATA (in matrix Xmat_c)
% covmat=cov(Xmat_c');
xc = Xmat_c';
xc_mean = mean(xc,1);
for n_ch=1:N_Channels
    xc(:,n_ch) = xc(:,n_ch) - xc_mean(n_ch);  % Remove mean
end
covmat = (xc' * xc) / (N-1);
clear xc

%%
if ~complexmixing, covmat=real(covmat); end
[Ec, Dc] = eig(covmat);
[dummy,order] = sort(diag(-Dc));
Ec = Ec(:,order(1:pcadim));
d = diag(Dc);
dsqrtinv = real(d.^(-0.5));
Dsqrtinv = diag(dsqrtinv(order(1:pcadim)));
whiteningmatrix = Dsqrtinv*Ec';
dewhiteningmatrix=Ec*inv(Dsqrtinv);
%THIS MATRIX IS THE MAIN INPUT TO THE ITERATIVE ALGORITHMS:
Zmat_c=whiteningmatrix*Xmat_c;

clear Xmat_c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FASTICA ESTIMATION

disp('Launching complex-valued FastICA')

%define constant in objective function. 1 or 0.1 seems to be the same
lambda=1;

%initial point, make it imaginary and unitary
W=randn(components,pcadim)...
    +complexmixing*sqrt(-1)*randn(components,pcadim);
W=inv(sqrtm(W*W'))*W;

%control when stopped
iter=0;
converged=false;

numofbits_maxiter = num2str(length(num2str(maxiter)));
fprintf('Iteration')
while ~converged
    %increase counter
    iter=iter+1;
    %Show progress to user
    eval(['fprintf('' %0',numofbits_maxiter,'d (max %0',numofbits_maxiter,'d)'',iter,maxiter);'])
    %store old value
    W_old=W;
    %Compute outputs, note lack of conjugate
    Y=W*Zmat_c;
    %Computing nonlinearities
    Y2=abs(Y).^2;
    gY2 = 1./(lambda + Y2);
    gderY2 = -gY2.^2;
    %Fixed-point iteration
    W=((Y.*gY2)*Zmat_c'/N - diag(mean(gY2 + Y2.*gderY2,2)) * W);
    %In case we want to restrict W to be real-valued, do it here:
    if ~complexmixing
        W=real(W);
    end
    %Make unitary
    if rcond(W*W')>1e-3
        W=sqrtm(inv(W*W'))*W;
    else
        W=sqrtm(pinv(W*W'))*W;
    end

    %check if converged
    convcriterion=1-trace(abs(W*W_old'))/components;
    if iter==maxiter | convcriterion<1e-7
        converged=1;
    else
        for nn=1:(length(num2str(maxiter))*2+8); fprintf('\b'); end        
    end
end
clear Y Y2 gY2 gderY2
fprintf('\n')
%Compute mixing matrix (in whitened space)
A=W';
%Compute source signal estimates
S=W*Zmat_c;
clear Zmat_c

%% SORT COMPONENTS AND TRANSFORM TO ORIGINAL SPACE
%compute objective function for each component
objective=[];
for i=1:size(W,1)
    objective(i)=-mean(log(lambda+(abs(S(i,:)).^2)));
end
%sort components using the objective
[dummy,componentorder]=sort(objective,'descend');
W=W(componentorder,:);
A=A(:,componentorder);
S=S(componentorder,:);
objective=objective(componentorder);

%Compute mixing and demixing matrix in original channel space
%Spatial filters (demixing matrix)
W_orig=W*whiteningmatrix;
%Spatial patterns (mixing matrix)
A_orig=dewhiteningmatrix*A;
%Independent components, but note that these are STFT's of the sources
%Adjust the dimension to (ic_chan,freq x time x trial)
for n_ic=1:components
    S_orig(n_ic,:,:,:) = reshape(S(n_ic,:),N_Freq,N_Time,N_Trials);
end
