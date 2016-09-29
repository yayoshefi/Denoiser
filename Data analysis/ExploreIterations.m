function []=ExploreIterations(iterations,Image,Noise)
global Analysis

if nargin==0;   I=Analysis.iterations;
else            I=iterations;
end

M=length(I);        step=2;%Analysis.DebuggerIter;
rows=3;

figure;     subplot (rows, ceil((M/step)/rows), 1)
imshow(I(1).Lhat,[]);      title('no context')
for t= step:step:length (I)
    subplot(rows, ceil((M/step)/rows), floor(t/step)+1)
    imshow(I(t).Lhat,[]);  title( sprintf('iteration %i',t) )
    
    cnt  = [I(t+1-step:t).changes];
    score= [I(t+1-step:t).score  ];
    l_0=I(t).l_0;
    sparsity= (l_0/(Analysis.K^2));
    
    if nargin >= 2
            [D_noised]=removenoise(double(Image),Noise,I(t).Lhat(:)');
            I(t).Psnr=psnr(D_noised,double(Image),255);
    end
    
%     fprintf('\niter %i: l_0 norm= %3.0f (%0.3f)\n',t,l_0,sparsity);
%     disp( strrep(['pixels changed:' sprintf(' %i,'  ,cnt(:) )], ',', '    '))
%     disp( strrep(['score  changed:' sprintf(' %3.0f,',score) ], ',', '    '))
end
colormap (Analysis.ColorMap);
disp('')
disp ( struct2table( rmfield(I,{'Lhat','CoOc'}) ) )
end

% step*floor(t/step)
% imshow( reshape( uint16(Lhat) , Analysis.LabelsSize),[] ); colormap (Analysis.ColorMap);