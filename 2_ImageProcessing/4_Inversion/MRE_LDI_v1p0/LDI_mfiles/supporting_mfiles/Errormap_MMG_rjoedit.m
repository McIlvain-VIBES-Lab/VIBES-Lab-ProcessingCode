function [U,ErrorMap]=Errormap_MMG_rjoedit(P);
% this code is taken from Strain_SNR_from_phase_MMG_rjoedit and is used to
% calculate an error map for uncertainty estimate
s=size(P);
U=zeros(s(1),s(2),s(3),3);
nph =s(4);  %number of phases
ErrorMap=zeros(s(1),s(2),s(3),s(5));
nph=s(4);
for jj=1:3
    FFTU=fft(P(:,:,:,:,jj),[],4); % take harmonic in time domain, same as isolate_harmonic
    U(:,:,:,jj)=2/nph*FFTU(:,:,:,2);
    t=(0:nph-1)*2*pi/nph;
    % Compute the disagreement of the harmonic estimate with the phase data as
    % an error estimate
    Ustack4D=P(:,:,:,:,jj);
    Ustack4D=Ustack4D-repmat(mean(Ustack4D,4),[1 1 1 nph]); % Remove DC component
    Uharm=zeros(s(1:4));
    for kk=1:nph
        Uharm(:,:,:,kk)=real(U(:,:,:,jj)*exp(1i*t(kk)));
    end
    ErrorMap(:,:,:,jj)=nph/(nph-3)*std(Ustack4D-Uharm,0,4);   % Apply correction because this error metric misses the noise in the DC and 1st harmonic.
end
Ur=real(U);
Ui=imag(U);
end