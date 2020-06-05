clear all
clc
Sigs = dir(['...\' '/' '*.mat']);
NumSigs = size(Sigs,1);
for i=1:NumSigs
    load(['....\',Sigs(i).name])
    zz=0;
    lim=1;
    final1=[];
    for qiu=1:17
        Samrate = 1000
        LEw = 1; 
        HEw = 75; 
        filterOrder = 2; 
        [b, a] = butter(filterOrder, [LEw HEw]/(Samrate/2)); 
        prdada = filtfilt(b, a, EEG.data(:,qiu)); 
        Soa =  1000;
        Sob = 200; 
        [N,D] = rat(Sob/Soa);
        Check = [Sob/Soa, N/D]; 
        S_200 = resample(prdada, N, D);
        [C,L] = wavedec(S_200, 5, 'db2');
        D1 = wrcoef('d',C,L,'db2',1);
        D2 = wrcoef('d',C,L,'db2',2);
        D3 = wrcoef('d',C,L,'db2',3);
        D4 = wrcoef('d',C,L,'db2',4);
        D5 = wrcoef('d',C,L,'db2',5);
        A5 = wrcoef('a',C,L,'db2',5);
        for j=1:885  
            segment_D2=D2(320*j-319:320*j);
            segment_D3=D3(320*j-319:320*j);
            segment_D4=D4(320*j-319:320*j);
            segment_D5=D5(320*j-319:320*j);
            segment_A5=A5(320*j-319:320*j);
            mu1=mean(segment_D2');
            mu2=mean(segment_D3');
            mu3=mean(segment_D4');
            mu4=mean(segment_D5');
            mu5=mean(segment_A5');
            Ham1= (std(fft(segment_D2))).*sum(([1:size(segment_D2,1)]'-mu1).*fft(segment_D2));
            Ham2= (std(fft(segment_D3))).*sum(([1:size(segment_D3,1)]'-mu2).*fft(segment_D3));
            Ham3= (std(fft(segment_D4))).*sum(([1:size(segment_D4,1)]'-mu3).*fft(segment_D4));
            Ham4= (std(fft(segment_D5))).*sum(([1:size(segment_D5,1)]'-mu4).*fft(segment_D5));
            Ham5= (std(fft(segment_A5))).*sum(([1:size(segment_A5,1)]'-mu5).*fft(segment_A5));
            En1=sum(fft(segment_D2).^2);
            En2=sum(fft(segment_D3).^2);
            En3=sum(fft(segment_D4).^2);
            En4=sum(fft(segment_D5).^2);
            En5=sum(fft(segment_A5).^2);
            HH1=sum(fft(segment_D2)./(1+abs([1:size(segment_D2,1)]')));
            HH2=sum(segment_D3./(1+abs([1:size(segment_D3,1)]')));
            HH3=sum(segment_D4./(1+abs([1:size(segment_D4,1)]')));
            HH4=sum(segment_D5./(1+abs([1:size(segment_D5,1)]')));
            HH5=sum(segment_A5./(1+abs([1:size(segment_A5,1)]')));
            Ent1=-sum(segment_D2.*log(segment_D2));
            Ent2=-sum(segment_D3.*log(segment_D3));
            Ent3=-sum(segment_D4.*log(segment_D4));
            Ent4=-sum(segment_D5.*log(segment_D5));
            Ent5=-sum(segment_A5.*log(segment_A5));  
            %     LOO(qiu,j,1:20)=real(log([Ham1;Ham2;Ham3;Ham4;Ham5;En1;En2;En3;En4;En5;HH1;HH2;HH3;HH4;HH5;Ent1;Ent2;Ent3;Ent4;Ent5]));
            LOO(qiu,j,1:5)=real(log([En1;En2;En3;En4;En5]));
            %     [Ham1;Ham2;Ham3;Ham4;Ham5;En1;En2;En3;En4;En5;HH1;HH2;HH3;HH4;HH5;Ent1;Ent2;Ent3;Ent4;Ent5];
        end
    end
    for p1=1:17
        for p2=1:5
        M1=reshape(LOO(p1,:,p2),[1,885]);
        M2=movingmean(M1',25,[],1);
        LOO(p1,1:885,p2)=M2;
        end
    end
    save(['...\',Sigs(i).name],'LOO')
    clearvars -except Sigs NumSigs
end


