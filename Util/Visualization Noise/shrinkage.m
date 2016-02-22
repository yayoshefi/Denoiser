function SV=shrinkage(SV, thresh, sigma, wsize, method)

switch method,
    case 'HTHRESH',  %Hard thresholding
        sind=find(abs(SV)<(thresh));
        SV(sind)=0;  %shrinkage by hard thresholding

    case 'NULL_LAST_COEFF',   % nulling the last coefficients
        SVS=SV.*SV;
        SVS=flipud(cumsum(flipud(SVS)));
        indmap=(SVS>thresh*(sigma*wsize)^2);
        SV=SV.*indmap;
        %imagesc(SV); pause;


    case 'NULL_SMALLEST_COEFF',   %shrinkage by nulling the smallest coefficients
        SVS=SV.*SV;
        for i=1:size(SV,2),
            [B,sind]=sort(SVS(:,i));
            S=SV(sind,i);
            A=cumsum(B);
            ind=find(A<thresh*(sigma*wsize)^2);
            S(ind)=0;
            SV(sind,i)=S;
        end;
        % imagesc(abs(SV)); pause;

    case 'WIENNER',  %not implelented

        %     D=diag(S);
        %     nD=D.*1./(1+1000./(D));
        %     %plot(nD);
        %     %plot(D,'k'); hold on; plot(nD,'r'); hold off;
        %     SV=diag(nD)*V';


end;