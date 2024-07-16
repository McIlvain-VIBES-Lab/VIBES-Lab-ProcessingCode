ti = 7;
threshold = 1;
waves = 2;
directions = 2;

Amp = [];
Adirs = [];
b = [];
clear Atot
if threshold == 1
    SF_thresh = 50;
    FA_thresh = .5;
else
    SF_thresh = 50;
    FA_thresh = 0;
end
FAv = FAtractlistAP(:,ti)>FA_thresh;
clear AtotS AtotF
for dirs = 1:length(Mutotlist(1,1,:))
    
    clear Mu
    Mu = Mutotlist(logical(FAv),ti,dirs);
    AwvS = []; AwvF = []; AMP = [];    
    for wv = 1:waves
        clear Sv Fv AS AF AmpS AmpF ArepS ArepF  
        
        AS(1:length(Mu),1) = 1;
        AS(1:length(Mu),2) = cos(Thlist(logical(FAv),wv,ti,dirs)).^2;
        AS(1:length(Mu),3) = 0;
        
        AF(1:length(Mu),1) = 1;
        AF(1:length(Mu),2) = cos(2*Thlist(logical(FAv),wv,ti,dirs)).^2;
        AF(1:length(Mu),3) = sin(2*Thlist(logical(FAv),wv,ti,dirs)).^2;
        
        Amp = Amptotlist(logical(FAv),ti,wv,dirs);
        
        ArepS = AS.*repmat(squeeze(Uperctot(logical(FAv),ti,((2*wv)-1),dirs)),[1 3])/100;
        ArepF = AF.*repmat(squeeze(Uperctot(logical(FAv),ti,(2*wv),dirs)),[1 3])/100;
        AwvS = cat(3,AwvS,ArepS); 
        AwvF = cat(3,AwvF,ArepF);
            
        AMP = cat(3,Amp,AMP);
        
    end
    if waves == 1
        Atot(1:length(Adirs),:,:,dirs) = Adirs;
    else
        AtotS(1:length(AwvS(:,1,1)),:,:,dirs) = AwvS.*(repmat(AMP,[1 3 1]));
        AtotF(1:length(AwvF(:,1,1)),:,:,dirs) = AwvF.*(repmat(AMP,[1 3 1]));
    end
    %sum(Mu>0)
    b(dirs,:) = Mu';
    %sum(b>0,2)
end

AtotS = permute(AtotS,[4 2 1 3]);
AtotF = permute(AtotF,[4 2 1 3]);

Aap = AtotS(1,:,:,1)+AtotS(1,:,:,2)+AtotF(1,:,:,1)+AtotF(1,:,:,2);
Alr = AtotS(2,:,:,1)+AtotS(2,:,:,2)+AtotF(2,:,:,1)+AtotF(2,:,:,2);
A = cat(1,Aap,Alr);


for ii = 1:length(A(1,1,:))
    x = squeeze(A(:,:,ii))\b(:,ii);
    bnew(:,ii) = A(:,:,ii)*x;
    x_soln(ii,:) = x';
end
